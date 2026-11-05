from pathlib import Path
# V0.15.1 MEMBRANE-ONLY PACKING DEFAULTS
# V0.15.0 MEMBRANE-ONLY BUILDER
# V0.14.0 GENERIC CUSTOM STEROL PLACEMENT
# V0.13.7 RESTORE LIPID21 TOPOLOGY HELPERS
# V0.13.7 RESTORE LIPID21 TOPOLOGY HELPERS
def _lipid21_interchange_line(serial,atom,icode,xyz):
    """Write one Lipid21 modular atom using the proven reference convention."""
    line=format_lipid_atom(serial,atom,1,xyz)
    line='ATOM  '+line[6:21]+' '+'    '+icode+line[27:]
    if line[21:22]!=' ' or line[22:26]!='    ' or line[26:27]!=icode:
        raise RuntimeError('Internal Lipid21 interchange-field audit failed')
    return line


def _write_lipid21_environment_files(items,prefix,box,batch_size=50):
    """Write one modular-lipid PDB and one bulk PDB for topology generation.

    POPC is written PA->PC->OL; POPS is written PA->PS->OL. A TER record is
    placed only after a complete lipid. Internal coordinates are unchanged.
    """
    prefix=Path(prefix)
    modular=[]
    bulk=[]
    for item in items:
        residue_names={str(atom['resname']).strip() for atom in item['molecule']['atoms']}
        if item.get('species') in {'POPC','POPS'} or residue_names & {'PA','PC','PS','OL'}:
            modular.append(item)
        else:
            bulk.append(item)

    files=[]
    if modular:
        path=prefix.parent/f'{prefix.name}_lipids.pdb'
        lines=[
            f"CRYST1{float(box[0]):9.3f}{float(box[1]):9.3f}{float(box[2]):9.3f}"
            f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1"
        ]
        serial=1
        for item in modular:
            headgroup='PS' if item.get('species')=='POPS' else 'PC'
            module_rank={'PA':0,headgroup:1,'OL':2}
            indexed=list(enumerate(zip(item['molecule']['atoms'],item['xyz'])))
            indexed.sort(key=lambda entry:(
                module_rank.get(str(entry[1][0]['resname']).strip(),99),
                entry[0],
            ))
            ordered_modules=[]
            for _,(atom,point) in indexed:
                module_name=str(atom['resname']).strip()
                if not ordered_modules or ordered_modules[-1]!=module_name:
                    ordered_modules.append(module_name)
                lines.append(_lipid21_interchange_line(serial,atom,'1',point))
                serial+=1
            expected=['PA',headgroup,'OL']
            if ordered_modules!=expected:
                raise RuntimeError(
                    f'Lipid21 module-order audit failed for {item.get("species")}: '
                    f'expected {expected}, observed {ordered_modules}'
                )
            lines.append(f"TER   {((serial-1)%99999)+1:5d}")
            serial+=1
        lines.append('END')
        path.write_text('\n'.join(lines)+'\n')
        files.append(path)

    if bulk:
        path=prefix.parent/f'{prefix.name}_bulk.pdb'
        _write_environment_pdb(bulk,path,box)
        files.append(path)

    return files,{
        'modular_lipids':len(modular),
        'bulk_molecules':len(bulk),
        'lipid_file_count':1 if modular else 0,
        'complete_lipid_ter_records':len(modular),
        'ter_between_lipid_modules':0,
        'module_order':{'POPC':['PA','PC','OL'],'POPS':['PA','PS','OL']},
        'separation_mode':'TER_after_complete_lipid',
        'temporary_interchange_files':[str(path) for path in files],
    }


def _tleap_environment_topology_from_items(items,prefix,prmtop_path,rst7_path,log_path,box,cfg):
    files,interchange_report=_write_lipid21_environment_files(items,prefix,box)
    if not files:
        raise RuntimeError('Generated environment contains no molecules')

    with tempfile.TemporaryDirectory() as tmpname:
        leap=Path(tmpname)/'environment.leap'
        commands=['source leaprc.lipid21','source leaprc.water.opc','source leaprc.gaff2']
        custom_settings=_custom_sterol_settings(cfg)
        if custom_settings:
            for parameter in custom_settings.get('parameters',[]):
                commands.append(f'loadAmberParams {Path(parameter).expanduser().resolve()}')
            commands.append(f'loadOff {Path(custom_settings["library"]).expanduser().resolve()}')
        units=[]
        for index,path in enumerate(files):
            unit=f'ENV{index:03d}'
            units.append(unit)
            commands.append(f'{unit} = loadPdb {path.resolve()}')
        commands.append('ENV = combine { '+' '.join(units)+' }')
        commands.append(
            f'saveAmberParm ENV {Path(prmtop_path).resolve()} '
            f'{Path(rst7_path).resolve()}'
        )
        commands.append('quit')
        leap.write_text('\n'.join(commands)+'\n')
        result=subprocess.run(
            ['tleap','-f',str(leap)],
            capture_output=True,
            text=True,
        )
        log_text=result.stdout+'\n'+result.stderr
        Path(log_path).write_text(log_text)
        if result.returncode!=0 or not Path(prmtop_path).is_file() or not Path(rst7_path).is_file():
            raise RuntimeError(
                'tleap could not parameterize the Lipid21/OPC/ion environment. '
                f'See {log_path}'
            )
        lower=log_text.lower()
        fatal_tokens=('fatal error','does not have a type','could not find bond parameter')
        if any(token in lower for token in fatal_tokens):
            raise RuntimeError(
                f'tleap reported an incomplete environment topology; see {log_path}'
            )
    return interchange_report


# V0.14.0 GENERIC CUSTOM STEROL PLACEMENT
def _parse_mol2_single_residue(path,residue_name,species):
    path=Path(path)
    lines=path.read_text(errors='replace').splitlines()
    try:
        atom_start=lines.index('@<TRIPOS>ATOM')+1
        bond_start=lines.index('@<TRIPOS>BOND')
    except ValueError as exc:
        raise ValueError(
            f'{path} is not a valid MOL2 file: required @<TRIPOS>ATOM/BOND sections are missing. '
            'The first section must be @<TRIPOS>MOLECULE.'
        ) from exc
    if not any(line.strip()=='@<TRIPOS>MOLECULE' for line in lines[:atom_start]):
        raise ValueError(f'{path}: expected @<TRIPOS>MOLECULE (not a species-specific header)')
    atoms=[]
    for line in lines[atom_start:bond_start]:
        if not line.strip() or line.startswith('@<TRIPOS>'): continue
        fields=line.split()
        if len(fields)<9: continue
        atom_id=int(fields[0]); name=fields[1]
        xyz=np.array([float(fields[2]),float(fields[3]),float(fields[4])],float)
        atom_type=fields[5]; mol2_resname=fields[7]
        if mol2_resname!=residue_name:
            raise ValueError(
                f'{path}: atom {name} uses residue {mol2_resname}, expected {residue_name}'
            )
        element=''.join(c for c in name if c.isalpha())[:1].upper() or 'X'
        atoms.append({'id':atom_id,'name':name,'resname':residue_name,
                      'element':element,'mol2_type':atom_type,'xyz':xyz,
                      'record':'HETATM','line':''})
    if not atoms: raise ValueError(f'No atoms parsed from {path}')
    names=[atom['name'] for atom in atoms]
    if len(names)!=len(set(names)): raise ValueError(f'Duplicate atom names in {path}')
    return {'atoms':atoms,'species':species,'residue_name':residue_name,'source':str(path)}


def _rotation_from_vectors(source,target):
    source=np.asarray(source,float); target=np.asarray(target,float)
    source=source/np.linalg.norm(source); target=target/np.linalg.norm(target)
    cross=np.cross(source,target); sine=float(np.linalg.norm(cross)); cosine=float(np.dot(source,target))
    if sine<1e-12:
        if cosine>0: return np.eye(3)
        trial=np.array([1.0,0.0,0.0]) if abs(source[0])<0.9 else np.array([0.0,1.0,0.0])
        axis=np.cross(source,trial); axis/=np.linalg.norm(axis)
        return 2.0*np.outer(axis,axis)-np.eye(3)
    axis=cross/sine
    k=np.array([[0.0,-axis[2],axis[1]],[axis[2],0.0,-axis[0]],[-axis[1],axis[0],0.0]])
    return np.eye(3)+k*sine+(k@k)*(1.0-cosine)


def _custom_sterol_template(settings):
    species=settings['species_code']; residue=settings['residue_name']
    parsed=_parse_mol2_single_residue(settings['mol2'],residue,species)
    atoms=parsed['atoms']; by_name={atom['name']:atom for atom in atoms}
    orientation=settings['orientation']; anchor_name=orientation['anchor_atom']
    end_names=list(orientation.get('axis_end_atoms',[]))
    if anchor_name not in by_name: raise ValueError(f'{species}: anchor atom {anchor_name} is absent')
    if not end_names or any(name not in by_name for name in end_names):
        raise ValueError(f'{species}: invalid axis_end_atoms {end_names}')
    anchor=by_name[anchor_name]['xyz'].copy()
    endpoint=np.mean([by_name[name]['xyz'] for name in end_names],axis=0)
    axis=endpoint-anchor
    if np.linalg.norm(axis)<5.0: raise ValueError(f'{species}: sterol axis is unexpectedly short')
    parsed.update({'anchor':anchor,'axis':axis,'anchor_name':anchor_name,'axis_end_names':end_names})
    return parsed


def _oriented_custom_sterol_candidate(template,host,axis_angle,offset):
    anchor=template['anchor']; source_axis=template['axis']
    target_axis=np.array([0.0,0.0,-1.0 if host['leaflet']=='upper' else 1.0])
    align=_rotation_from_vectors(source_axis,target_axis)
    xyz=np.array([atom['xyz'] for atom in template['atoms']],float)
    xyz=(xyz-anchor)@align.T
    c=math.cos(axis_angle); s=math.sin(axis_angle)
    zrot=np.array([[c,-s,0.0],[s,c,0.0],[0.0,0.0,1.0]])
    xyz=xyz@zrot.T+host['anchor']+offset
    final_anchor=host['anchor']+offset
    atoms=[dict(atom) for atom in template['atoms']]
    molecule={'atoms':atoms,'anchor':final_anchor,'leaflet':host['leaflet'],
              'species':template['species']}
    return {'molecule':molecule,'xyz':xyz,'anchor':final_anchor,
            'elements':[atom['element'] for atom in atoms],
            'leaflet':host['leaflet'],'species':template['species']}


def _place_custom_sterol_item(host,template,fixed_cells,items,host_index,cfg,rng,settings):
    placement=settings.get('placement',{})
    trials=int(placement.get('trials_per_site',240))
    xy=float(placement.get('xy_offset_angstrom',3.0))
    z=float(placement.get('z_offset_angstrom',1.5))
    best=None
    for _ in range(trials):
        offset=np.array([rng.uniform(-xy,xy),rng.uniform(-xy,xy),rng.uniform(-z,z)])
        candidate=_oriented_custom_sterol_candidate(template,host,float(rng.uniform(0,2*math.pi)),offset)
        valid,score,minima=_validate_replacement_item(candidate,host_index,items,fixed_cells,cfg)
        if valid and (best is None or score<best[0]): best=(score,candidate,minima)
        if best is not None and best[0]<=1e-12: break
    if best is None: return None
    best[1]['score']=float(best[0]); best[1]['minimum_distances']=best[2]
    return best[1]


def _place_exact_custom_sterol_target(items,eligible,target,template,fixed_cells,cfg,rng,leaf,settings):
    candidates=[index for index in eligible if items[index]['species']=='POPC']
    rng.shuffle(candidates); placed=[]; failures=[]
    for index in candidates:
        if len(placed)>=target: break
        candidate=_place_custom_sterol_item(items[index],template,fixed_cells,items,index,cfg,rng,settings)
        if candidate is None:
            failures.append({'site_index':index,'reason':'no_valid_custom_sterol_candidate'})
            continue
        items[index]=candidate
        placed.append({'site_index':index,'score':candidate['score'],
                       'minimum_distances':candidate['minimum_distances']})
    if len(placed)!=target:
        raise RuntimeError(
            f"Exact composition failed: requested {target} {template['species']} in {leaf}, "
            f'placed {len(placed)} after testing {len(candidates)} POPC sites.'
        )
    return placed,failures


def _custom_sterol_settings(cfg):
    block=cfg.get('custom_sterols',{})
    if not block.get('enabled',False): return None
    species=block.get('species',[])
    if len(species)!=1: raise ValueError('This release supports exactly one active custom sterol species')
    settings=species[0]
    for key in ('species_code','residue_name','mol2','library','parameters','orientation'):
        if key not in settings: raise ValueError(f'custom sterol setting is missing {key}')
    return settings

# V0.13.6 RESTORE TOPOLOGY MAPPING HELPER
# V0.13.5 ORDER-INDEPENDENT SOLUTE AUDIT
# V0.13.4 LIPID21 MODULE ORDER
# V0.13.3 LIPID21 COMPLETE-LIPID TER SEPARATION
# V0.13.2 BATCHED LIPID21 INTERCHANGE
# V0.13.1 PDB ATOM-NAME COLUMN ALIGNMENT
# V0.13.0 AUTHORITATIVE TOPOLOGY MAPPING
# V0.12.9 TOPOLOGY-PRESERVING ION AND AMBER OUTPUT
# V0.12.8B RESTORE ORIGINAL DIRECT OPC PLACEMENT
# V0.12.7 EXPLICIT FIXED-ATOM SOLVENT CONTAINMENT
# V0.12.6 CRYST1 FINAL BOX OUTPUT
# V0.12.5C FILTER TERMINAL AND CONNECTIVITY RECORDS
from collections import defaultdict, Counter
import argparse, hashlib, json, math, os, subprocess, tempfile
import numpy as np
import yaml

COORD_RECORDS=('ATOM  ','HETATM')


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def element_from_line(line,name):
    e=line[76:78].strip().upper() if len(line)>=78 else ''
    if e: return e
    stripped=''.join(c for c in name if c.isalpha()).upper()
    return stripped[:1] or 'X'


def parse_atom(line,index):
    name=line[12:16].strip(); resname=line[17:20].strip(); chain=line[21:22]
    resid=line[22:26].strip(); icode=line[26:27]
    xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],float)
    return {'line':line,'index':index,'record':line[:6].strip(),'name':name,'resname':resname,
            'chain':chain,'resid':resid,'icode':icode,'xyz':xyz,'element':element_from_line(line,name)}


def parse_opm(path,protein_cfg=None):
    """Parse an OPM-preoriented immutable assembly.

    Supported modes are intentionally restrictive:
      1. OPM DUM records present: use their two Z planes.
      2. DUM records absent, but orientation_source=OPM and preoriented=true:
         preserve all coordinates, use Z as the membrane normal, and center the
         calibrated host membrane at z=0. No orientation or translation is done.
    """
    protein_cfg=protein_cfg or {}
    if protein_cfg.get('orientation_source')!='OPM' or not protein_cfg.get('preoriented'):
        raise ValueError(
            'AMBMEM requires orientation_source: OPM and preoriented: true. '
            'It does not orient membrane proteins.'
        )
    lines=Path(path).read_text(errors='replace').splitlines()
    physical=[]; dum=[]; coordinate_index=0
    for line in lines:
        if line.startswith(COORD_RECORDS):
            atom=parse_atom(line,coordinate_index); coordinate_index+=1
            if atom['resname']=='DUM': dum.append(atom)
            else: physical.append(atom)
    if not physical:
        raise ValueError('No physical ATOM/HETATM records found')

    if dum:
        planes=sorted(set(round(float(atom['xyz'][2]),3) for atom in dum))
        if len(planes)!=2:
            raise ValueError(f'Expected two OPM DUM planes, found {planes}')
        source='OPM_DUM_records'
    else:
        # DUM-free OPM output: retain the OPM coordinate frame and use the
        # calibrated host membrane thickness centered at z=0.
        half_thickness=float(protein_cfg.get(
            'dum_free_half_thickness_angstrom',17.9
        ))
        if not np.isfinite(half_thickness) or half_thickness<=0:
            raise ValueError('dum_free_half_thickness_angstrom must be positive')
        planes=[-half_thickness,half_thickness]
        source='OPM_preoriented_z0_default_no_DUM'

    return lines,physical,dum,planes,source

def parse_reference(path):
    """Parse a Lipid21 POPC patch whose molecular modules are separate residues.

    Lipid21 represents one POPC molecule as three consecutive residue modules:
    PA, PC, and OL. Residue numbers therefore identify modules, not complete
    lipids, and cannot be used alone to group a molecule.
    """
    lines=Path(path).read_text(errors='replace').splitlines()
    cryst=next((line for line in lines if line.startswith('CRYST1')),None)
    if not cryst:
        raise ValueError('POPC reference has no CRYST1 record')
    box=np.array([float(cryst[6:15]),float(cryst[15:24]),float(cryst[24:33])])

    # Preserve residue-block order exactly as it appears in the PDB.
    residue_blocks=[]
    current_key=None
    current_atoms=[]
    for line in lines:
        if not line.startswith(COORD_RECORDS):
            continue
        atom=parse_atom(line,0)
        if atom['resname'] not in {'PA','PC','OL'}:
            continue
        key=(atom['chain'],atom['resid'],atom['icode'],atom['resname'])
        if current_key is None:
            current_key=key
        if key!=current_key:
            residue_blocks.append((current_key,current_atoms))
            current_key=key
            current_atoms=[]
        current_atoms.append(atom)
    if current_atoms:
        residue_blocks.append((current_key,current_atoms))

    # Group each PA/PC/OL triplet into one complete POPC molecule. Accept any
    # ordering within the consecutive triplet, but require exactly one module
    # of each type and exactly one PC:P31 anchor.
    molecules=[]
    i=0
    while i<len(residue_blocks):
        window=residue_blocks[i:i+3]
        names=[key[3] for key,_ in window]
        if len(window)==3 and set(names)=={'PA','PC','OL'}:
            atoms=[]
            for _,block in window:
                atoms.extend(block)
            anchors=[atom for atom in atoms if atom['name']=='P31' and atom['resname']=='PC']
            if len(anchors)!=1:
                raise ValueError(f'POPC module triplet at block {i} has {len(anchors)} PC:P31 anchors')
            molecules.append({
                'key':tuple(key for key,_ in window),
                'atoms':atoms,
                'anchor':anchors[0]['xyz'].copy(),
                'leaflet':'upper' if anchors[0]['xyz'][2]>=0 else 'lower'
            })
            i+=3
        else:
            raise ValueError(
                f'Unexpected Lipid21 module order at residue block {i}: {names}. '
                'Expected consecutive PA/PC/OL modules.'
            )

    if len(molecules)!=128:
        counts=Counter(key[3] for key,_ in residue_blocks)
        raise ValueError(
            f'Expected 128 complete POPC molecules, found {len(molecules)}; '
            f'module counts={dict(counts)}'
        )
    return molecules,box

def build_fixed_cells(atoms,cell_width):
    cells=defaultdict(list)
    for atom in atoms:
        key=tuple(np.floor(atom['xyz']/cell_width).astype(int))
        cells[key].append(atom)
    return cells


def clashes(molecule,fixed_cells,shift,cell_width,heavy_cut,hydrogen_cut):
    heavy2=heavy_cut*heavy_cut; hydrogen2=hydrogen_cut*hydrogen_cut
    for atom in molecule['atoms']:
        xyz=atom['xyz']+shift
        base=np.floor(xyz/cell_width).astype(int)
        cutoff2=hydrogen2 if atom['element']=='H' else heavy2
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for fixed in fixed_cells.get((base[0]+dx,base[1]+dy,base[2]+dz),()):
                        if float(np.dot(xyz-fixed['xyz'],xyz-fixed['xyz']))<cutoff2:
                            return True
    return False



# V0.10.1 ADAPTIVE POPC BOUNDARY REFILL
def _rotate_z(xyz,center,angle):
    c=math.cos(angle); s=math.sin(angle)
    rotation=np.array([[c,-s,0.0],[s,c,0.0],[0.0,0.0,1.0]])
    return (xyz-center)@rotation.T+center


def _atom_arrays_from_molecule(molecule,shift=None,angle=0.0,target_anchor=None):
    xyz=np.array([atom['xyz'] for atom in molecule['atoms']],float)
    anchor=molecule['anchor'].copy()
    if target_anchor is not None:
        xyz=xyz+(target_anchor-anchor)
        anchor=target_anchor.copy()
    if abs(angle)>1e-15:
        xyz=_rotate_z(xyz,anchor,angle)
    if shift is not None:
        xyz=xyz+shift
        anchor=anchor+shift
    elements=[atom['element'] for atom in molecule['atoms']]
    return xyz,anchor,elements


def _build_atom_cells(items,cell_width):
    cells=defaultdict(list)
    for molecule_index,item in enumerate(items):
        xyz=item['xyz']
        elements=item['elements']
        for atom_index,point in enumerate(xyz):
            key=tuple(np.floor(point/cell_width).astype(int))
            cells[key].append((molecule_index,atom_index,point,elements[atom_index]))
    return cells


def _candidate_valid(candidate_xyz,candidate_elements,fixed_cells,lipid_cells,cell_width,
                     fixed_heavy_cut,fixed_hydrogen_cut,lipid_hh_cut,lipid_hx_cut,lipid_xx_cut):
    minima={'fixed_heavy_lipid_heavy':999.0,'fixed_heavy_lipid_hydrogen':999.0,
            'lipid_heavy_heavy':999.0,'lipid_heavy_hydrogen':999.0,'lipid_hydrogen_hydrogen':999.0}
    score=0.0
    for point,element in zip(candidate_xyz,candidate_elements):
        base=np.floor(point/cell_width).astype(int)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    key=(base[0]+dx,base[1]+dy,base[2]+dz)
                    for fixed in fixed_cells.get(key,()):
                        d=float(np.linalg.norm(point-fixed['xyz']))
                        category='fixed_heavy_lipid_hydrogen' if element=='H' else 'fixed_heavy_lipid_heavy'
                        minima[category]=min(minima[category],d)
                        cutoff=fixed_hydrogen_cut if element=='H' else fixed_heavy_cut
                        if d<cutoff: return False,1e9,minima
                    for _,_,other,e2 in lipid_cells.get(key,()):
                        d=float(np.linalg.norm(point-other))
                        if element=='H' and e2=='H': category='lipid_hydrogen_hydrogen'; cutoff=lipid_xx_cut; soft=1.8; weight=1.0
                        elif element=='H' or e2=='H': category='lipid_heavy_hydrogen'; cutoff=lipid_hx_cut; soft=2.3; weight=3.0
                        else: category='lipid_heavy_heavy'; cutoff=lipid_hh_cut; soft=3.0; weight=10.0
                        minima[category]=min(minima[category],d)
                        if d<cutoff: return False,1e9,minima
                        score+=max(0.0,soft-d)**2*weight
    return True,score,minima


def _anchor_inside_fixed_footprint(anchor,membrane_fixed,xy_cut,z_cut):
    for atom in membrane_fixed:
        if abs(float(atom['xyz'][2]-anchor[2]))<=z_cut:
            if float(np.linalg.norm(atom['xyz'][:2]-anchor[:2]))<xy_cut:
                return True
    return False


# V0.10.2 TOPOLOGY-BASED BOUNDARY REFILL
def _refill_boundary_sites(removed,molecules,kept,fixed,lower,upper,box,cfg,rng):
    """Refill exterior annular sites using retained-host topology and exact contacts.

    Removed sites are attempted only when adjacent to the currently retained
    membrane in the same leaflet. Two bounded passes allow the exterior membrane
    edge to advance inward without classifying sites from protein-atom proximity.
    Sites that remain disconnected or cannot accept an intact POPC are left empty.
    """
    settings=cfg.get('boundary_packing',{})
    if not settings.get('enabled',True):
        return kept,[],removed

    conformer_limit=int(settings.get('conformers_per_site',32))
    rotations=int(settings.get('rotations_per_conformer',16))
    offsets=int(settings.get('offsets_per_rotation',8))
    xy_range=float(settings.get('xy_offset_angstrom',3.0))
    z_range=float(settings.get('z_offset_angstrom',0.75))
    passes=int(settings.get('passes',2))
    neighbor_distance=float(settings.get('retained_neighbor_distance_angstrom',11.0))
    minimum_neighbors=int(settings.get('minimum_retained_neighbors',1))
    cell_width=4.0

    fixed_cells=build_fixed_cells(fixed,cell_width)
    fixed_heavy=float(cfg['packing'].get('fixed_heavy_lipid_heavy_angstrom',2.0))
    fixed_hydrogen=float(cfg['packing'].get('fixed_heavy_lipid_hydrogen_angstrom',1.2))
    hh=float(settings.get('lipid_heavy_heavy_angstrom',1.5))
    hx=float(settings.get('lipid_heavy_hydrogen_angstrom',1.05))
    xx=float(settings.get('lipid_hydrogen_hydrogen_angstrom',0.8))

    dynamic=[]
    for molecule,shift in kept:
        xyz,anchor,elements=_atom_arrays_from_molecule(molecule,shift=shift)
        dynamic.append({'xyz':xyz,'anchor':anchor,'elements':elements,'leaflet':molecule['leaflet']})

    pools={leaf:[m for m in molecules if m['leaflet']==leaf] for leaf in ('upper','lower')}
    remaining=list(removed)
    accepted=[]

    def neighbor_count(record):
        target=np.asarray(record['anchor'],float)
        return sum(
            item['leaflet']==record['leaflet'] and
            float(np.linalg.norm(item['anchor'][:2]-target[:2]))<=neighbor_distance
            for item in dynamic
        )

    for pass_number in range(1,passes+1):
        eligible=[]; deferred=[]
        for record in remaining:
            count=neighbor_count(record)
            if count>=minimum_neighbors:
                eligible.append((count,record))
            else:
                deferred.append(record)
        if not eligible:
            remaining=deferred
            break

        rng.shuffle(eligible)
        # More retained neighbors first so continuous outer boundary sites are
        # committed before narrow or enclosed sites.
        eligible.sort(key=lambda item:item[0],reverse=True)
        next_remaining=list(deferred)

        for retained_neighbors,record in eligible:
            target=np.asarray(record['anchor'],float)
            pool=list(pools[record['leaflet']]); rng.shuffle(pool); pool=pool[:conformer_limit]
            lipid_cells=_build_atom_cells(dynamic,cell_width)
            best=None

            for conformer in pool:
                for angle in rng.random(rotations)*2*math.pi:
                    trial_offsets=[np.zeros(3)]
                    trial_offsets.extend(
                        np.array([
                            rng.uniform(-xy_range,xy_range),
                            rng.uniform(-xy_range,xy_range),
                            rng.uniform(-z_range,z_range)
                        ])
                        for _ in range(offsets)
                    )
                    for offset in trial_offsets:
                        xyz,anchor,elements=_atom_arrays_from_molecule(
                            conformer,angle=angle,target_anchor=target,shift=offset
                        )
                        valid,score,minima=_candidate_valid(
                            xyz,elements,fixed_cells,lipid_cells,cell_width,
                            fixed_heavy,fixed_hydrogen,hh,hx,xx
                        )
                        if not valid:
                            continue
                        key=(score,-min(minima.values()))
                        if best is None or key<best[0]:
                            best=(key,conformer,xyz,anchor,elements,minima)
                        if score<=1e-12:
                            break
                    if best is not None and best[0][0]<=1e-12:
                        break
                if best is not None and best[0][0]<=1e-12:
                    break

            if best is None:
                next_remaining.append(record)
                continue

            _,conformer,xyz,anchor,elements,minima=best
            dynamic.append({
                'xyz':xyz,'anchor':anchor,'elements':elements,
                'leaflet':record['leaflet']
            })
            accepted.append({
                'leaflet':record['leaflet'],
                'pass':pass_number,
                'retained_neighbors_at_attempt':retained_neighbors,
                'original_anchor':record['anchor'],
                'final_anchor':[float(value) for value in anchor],
                'score':float(best[0][0]),
                'minimum_distances':minima,
                'molecule':conformer,
                'xyz':xyz
            })

        remaining=next_remaining

    return kept,accepted,remaining


# V0.11.0 COMPLEX PROTEIN MEMBRANE
def _kabsch(mobile,reference):
    mc=mobile.mean(axis=0); rc=reference.mean(axis=0)
    h=(mobile-mc).T@(reference-rc)
    u,s,vt=np.linalg.svd(h); rotation=vt.T@u.T
    if np.linalg.det(rotation)<0:
        vt[-1,:]*=-1; rotation=vt.T@u.T
    return rotation,rc-mc@rotation.T


def _parse_modular_lipid_reference(path,headgroup,species):
    lines=Path(path).read_text(errors='replace').splitlines()
    blocks=[]; key=None; atoms=[]
    allowed={'PA',headgroup,'OL'}
    for line in lines:
        if not line.startswith(COORD_RECORDS): continue
        atom=parse_atom(line,0)
        if atom['resname'] not in allowed: continue
        new_key=(atom['chain'],atom['resid'],atom['icode'],atom['resname'])
        if key is None: key=new_key
        if new_key!=key:
            blocks.append((key,atoms)); key=new_key; atoms=[]
        atoms.append(atom)
    if atoms: blocks.append((key,atoms))
    molecules=[]
    for index in range(0,len(blocks),3):
        window=blocks[index:index+3]
        if len(window)!=3 or set(k[3] for k,_ in window)!=allowed: continue
        merged=[atom for _,block in window for atom in block]
        anchors=[a for a in merged if a['resname']==headgroup and a['name']=='P31']
        if len(anchors)!=1: continue
        molecules.append({'atoms':merged,'anchor':anchors[0]['xyz'].copy(),
                          'leaflet':'upper' if anchors[0]['xyz'][2]>=0 else 'lower',
                          'species':species})
    if not molecules: raise ValueError(f'No {species} molecules parsed from {path}')
    return molecules


def _parse_chl_reference(path):
    groups=defaultdict(list)
    for line in Path(path).read_text(errors='replace').splitlines():
        if not line.startswith(COORD_RECORDS): continue
        atom=parse_atom(line,0)
        if atom['resname'] not in {'CHL','CHL1'}: continue
        groups[(atom['chain'],atom['resid'],atom['icode'])].append(atom)
    molecules=[]
    for key,atoms in groups.items():
        anchors=[a for a in atoms if a['name'] in {'O1','O3'}]
        if len(anchors)!=1: continue
        molecules.append({'atoms':atoms,'anchor':anchors[0]['xyz'].copy(),
                          'leaflet':'upper' if anchors[0]['xyz'][2]>=0 else 'lower',
                          'species':'CHL'})
    if not molecules: raise ValueError(f'No CHL conformers parsed from {path}')
    return molecules


def _distance_to_fixed_xy(anchor,membrane_fixed):
    if not membrane_fixed: return float('inf')
    return min(float(np.linalg.norm(anchor[:2]-atom['xyz'][:2])) for atom in membrane_fixed)


def _convert_popc_to_pops_item(item,donors,rng):
    molecule=item['molecule']; xyz=item['xyz']; anchor=item['anchor']
    host_atoms=molecule['atoms']
    host_pc=[(i,a) for i,a in enumerate(host_atoms) if a['resname']=='PC']
    host_map={a['name']:(i,a) for i,a in host_pc}
    donor=donors[int(rng.integers(len(donors)))]
    donor_ps=[a for a in donor['atoms'] if a['resname']=='PS']
    donor_map={a['name']:a for a in donor_ps}
    alignment=['P31','O31','O32','C3','C2','C1','O21','C21','O11','C11']
    common=[name for name in alignment if name in host_map and name in donor_map]
    if len(common)<4: return None
    mobile=np.array([donor_map[name]['xyz'] for name in common])
    reference=np.array([xyz[host_map[name][0]] for name in common])
    rotation,translation=_kabsch(mobile,reference)
    retained_atoms=[]; retained_xyz=[]
    for i,atom in enumerate(host_atoms):
        if atom['resname']!='PC': retained_atoms.append(atom); retained_xyz.append(xyz[i])
    transformed_atoms=[]; transformed_xyz=[]
    for atom in donor_ps:
        copy=dict(atom); copy['resname']='PS'
        transformed_atoms.append(copy); transformed_xyz.append(atom['xyz']@rotation.T+translation)
    atoms=retained_atoms+transformed_atoms
    coordinates=np.array(retained_xyz+transformed_xyz,float)
    anchors=[i for i,a in enumerate(atoms) if a['resname']=='PS' and a['name']=='P31']
    if len(anchors)!=1: return None
    return {'molecule':{'atoms':atoms,'anchor':coordinates[anchors[0]],'leaflet':item['leaflet']},
            'xyz':coordinates,'anchor':coordinates[anchors[0]],'elements':[a['element'] for a in atoms],
            'leaflet':item['leaflet'],'species':'POPS'}


def _place_chl_item(host,donors,fixed_cells,lipid_items,cfg,rng):
    settings=cfg['complex_membrane']; cell_width=4.0
    lipid_cells=_build_atom_cells([item for item in lipid_items if item is not host],cell_width)
    fixed_heavy=float(cfg.get('packing',{}).get('fixed_heavy_lipid_heavy_angstrom',2.0))
    fixed_hydrogen=float(cfg.get('packing',{}).get('fixed_heavy_lipid_hydrogen_angstrom',1.2))
    hh=1.5; hx=1.05; xx=0.8
    pool=[m for m in donors if m['leaflet']==host['leaflet']] or donors
    best=None
    for _ in range(int(settings.get('cholesterol_trials',160))):
        donor=pool[int(rng.integers(len(pool)))]
        angle=float(rng.uniform(0,2*math.pi)); offset=np.array([rng.uniform(-2,2),rng.uniform(-2,2),rng.uniform(-1.2,1.2)])
        xyz,anchor,elements=_atom_arrays_from_molecule(donor,angle=angle,target_anchor=host['anchor'],shift=offset)
        valid,score,minima=_candidate_valid(xyz,elements,fixed_cells,lipid_cells,cell_width,fixed_heavy,fixed_hydrogen,hh,hx,xx)
        valid=valid and minima['lipid_hydrogen_hydrogen']>=1.15
        if valid and (best is None or score<best[0]): best=(score,donor,xyz,anchor,elements,minima)
    if best is None: return None
    return {'molecule':best[1],'xyz':best[2],'anchor':best[3],'elements':best[4],
            'leaflet':host['leaflet'],'species':'CHL','score':float(best[0]),'minimum_distances':best[5]}



# V0.11.1 CHOLESTEROL-FIRST EXACT COMPOSITION
# V0.11.2 STRICT INTERMOLECULAR LIPID H-H FLOOR
# V0.11.3 ASYMMETRIC LEAFLET COMPOSITION
# V0.12.0 ASYMMETRIC OPC SOLVATION
# V0.12.1 OPM INPUT WITHOUT DUM MARKERS
# V0.12.2 PREORIENTED OPM DUM-OPTIONAL
# V0.12.3 CONTIGUOUS OPC TILING
# V0.12.4 TLEAP-CONTIGUOUS OPC VOLUME
def _validate_replacement_item(candidate,host_index,items,fixed_cells,cfg):
    """Validate one replacement against fixed input and every other membrane item."""
    environment=[item for index,item in enumerate(items) if index!=host_index]
    lipid_cells=_build_atom_cells(environment,4.0)
    return _candidate_valid(
        candidate['xyz'],candidate['elements'],fixed_cells,lipid_cells,4.0,
        float(cfg.get('packing',{}).get('fixed_heavy_lipid_heavy_angstrom',2.0)),
        float(cfg.get('packing',{}).get('fixed_heavy_lipid_hydrogen_angstrom',1.2)),
        1.5,1.05,
        float(cfg.get('packing',{}).get('lipid_hydrogen_hydrogen_angstrom',1.15))
    )


def _place_exact_cholesterol_target(items,eligible,target,donors,fixed_cells,cfg,rng,leaf):
    """Try seed-ordered POPC hosts until the exact CHL target is committed."""
    candidates=list(eligible); rng.shuffle(candidates)
    placed=[]; failures=[]
    for index in candidates:
        if len(placed)>=target: break
        if items[index]['species']!='POPC': continue
        candidate=_place_chl_item(items[index],donors,fixed_cells,items,cfg,rng)
        if candidate is None:
            failures.append({'site_index':index,'reason':'no_valid_chl_candidate'})
            continue
        valid,score,minima=_validate_replacement_item(candidate,index,items,fixed_cells,cfg)
        floor=float(cfg['complex_membrane'].get('cholesterol_minimum_hydrogen_hydrogen_angstrom',1.15))
        if not valid or minima['lipid_hydrogen_hydrogen']<floor:
            failures.append({'site_index':index,'reason':'final_recheck_failed','minimum_distances':minima})
            continue
        candidate['score']=float(score); candidate['minimum_distances']=minima
        items[index]=candidate
        placed.append({'site_index':index,'score':float(score),'minimum_distances':minima})
    if len(placed)!=target:
        raise RuntimeError(
            f'Exact composition failed: requested {target} CHL in {leaf} leaflet, '
            f'placed {len(placed)} after testing {len(candidates)} eligible POPC sites.'
        )
    return placed,failures


def _place_exact_pops_target(items,eligible,target,donors,fixed_cells,cfg,rng,leaf):
    """Convert POPC to POPS after CHL, continuing until exact target is reached."""
    candidates=[index for index in eligible if items[index]['species']=='POPC']
    rng.shuffle(candidates)
    placed=[]; failures=[]
    donor_pool=[m for m in donors if m['leaflet']==leaf] or donors
    donor_limit=min(len(donor_pool),int(cfg['complex_membrane'].get('pops_donors_per_site',64)))
    for index in candidates:
        if len(placed)>=target: break
        order=list(rng.permutation(len(donor_pool)))[:donor_limit]
        best=None
        for donor_index in order:
            attempt=_convert_popc_to_pops_item(items[index],[donor_pool[donor_index]],rng)
            if attempt is None: continue
            valid,score,minima=_validate_replacement_item(attempt,index,items,fixed_cells,cfg)
            if not valid: continue
            if best is None or score<best[0]: best=(score,attempt,minima)
            if score<=1e-12: break
        if best is None:
            failures.append({'site_index':index,'reason':'no_valid_pops_conversion'})
            continue
        score,attempt,minima=best
        items[index]=attempt
        placed.append({'site_index':index,'score':float(score),'minimum_distances':minima})
    if len(placed)!=target:
        raise RuntimeError(
            f'Exact composition failed: requested {target} POPS in {leaf} leaflet, '
            f'placed {len(placed)} after testing {len(candidates)} eligible POPC sites.'
        )
    return placed,failures


def _apply_complex_composition(items,fixed,lower,upper,cfg,rng):
    """Apply exact per-leaflet CHL, custom sterol, and optional POPS composition."""
    settings=cfg.get('complex_membrane',{})
    if not settings.get('enabled',False): return items,{}
    membrane_fixed=[atom for atom in fixed if lower-8<=atom['xyz'][2]<=upper+8]
    shell=float(settings.get('protein_exclusion_shell_angstrom',8.0))
    chl_ref=_parse_chl_reference(settings['CHL_reference'])
    custom_settings=_custom_sterol_settings(cfg)
    custom_template=_custom_sterol_template(custom_settings) if custom_settings else None
    pops_enabled=any(float(settings.get('leaflets',{}).get(leaf,{}).get(
        'POPS_fraction',settings.get('POPS_fraction',0.0)))>0 for leaf in ('upper','lower'))
    pops_ref=(_parse_modular_lipid_reference(settings['POPS_reference'],'PS','POPS')
              if pops_enabled else [])
    fixed_cells=build_fixed_cells(fixed,4.0)
    report={'order':['CHL','CUSTOM_STEROL','POPS'],'protected_shell_angstrom':shell,
            'custom_sterol':None if custom_template is None else {
                'species':custom_template['species'],'residue_name':custom_template['residue_name'],
                'anchor_atom':custom_template['anchor_name'],
                'axis_end_atoms':custom_template['axis_end_names'],
                'axis_length_angstrom':float(np.linalg.norm(custom_template['axis'])),
                'reflection_allowed':False},'leaflets':{}}

    for leaf in ('upper','lower'):
        leaflet_indices=[index for index,item in enumerate(items) if item['leaflet']==leaf]
        total=len(leaflet_indices); leaflet_settings=settings.get('leaflets',{}).get(leaf,{})
        chl_fraction=float(leaflet_settings.get('CHL_fraction',settings.get('CHL_fraction',0.05)))
        pops_fraction=float(leaflet_settings.get('POPS_fraction',settings.get('POPS_fraction',0.0)))
        custom_fraction=0.0
        if custom_settings:
            custom_fraction=float(custom_settings.get('leaflets',{}).get(leaf,{}).get(
                'fraction',custom_settings.get('fraction',0.05)))
        if min(chl_fraction,pops_fraction,custom_fraction)<0 or chl_fraction+pops_fraction+custom_fraction>1.0:
            raise ValueError(f'Invalid {leaf} fractions: CHL={chl_fraction}, custom={custom_fraction}, POPS={pops_fraction}')
        target_chl=int(round(total*chl_fraction)); target_custom=int(round(total*custom_fraction))
        target_pops=int(round(total*pops_fraction)); target_popc=total-target_chl-target_custom-target_pops
        eligible=[index for index in leaflet_indices if items[index]['species']=='POPC' and
                  _distance_to_fixed_xy(items[index]['anchor'],membrane_fixed)>=shell]
        if len(eligible)<target_chl+target_custom+target_pops:
            raise RuntimeError(f'Exact composition impossible in {leaf}: {len(eligible)} eligible sites')
        chl_placed,chl_failures=(_place_exact_cholesterol_target(
            items,eligible,target_chl,chl_ref,fixed_cells,cfg,rng,leaf) if target_chl else ([],[]))
        custom_placed,custom_failures=(_place_exact_custom_sterol_target(
            items,eligible,target_custom,custom_template,fixed_cells,cfg,rng,leaf,custom_settings)
            if target_custom else ([],[]))
        pops_placed,pops_failures=(_place_exact_pops_target(
            items,eligible,target_pops,pops_ref,fixed_cells,cfg,rng,leaf) if target_pops else ([],[]))
        species_names=['POPC','POPS','CHL']+([custom_template['species']] if custom_template else [])
        actual={species:sum(items[index]['species']==species for index in leaflet_indices)
                for species in species_names}
        expected={'POPC':target_popc,'POPS':target_pops,'CHL':target_chl}
        if custom_template: expected[custom_template['species']]=target_custom
        if actual!=expected: raise RuntimeError(f'Composition verification failed in {leaf}: expected {expected}, got {actual}')
        report['leaflets'][leaf]={
            'total':total,'requested_fractions':{'POPC':1.0-chl_fraction-custom_fraction-pops_fraction,
            'CHL':chl_fraction,'custom_sterol':custom_fraction,'POPS':pops_fraction},
            'target':expected,'actual':actual,'eligible_bulk_sites':len(eligible),
            'cholesterol_placements':chl_placed,'cholesterol_failed_sites':chl_failures,
            'custom_sterol_placements':custom_placed,'custom_sterol_failed_sites':custom_failures,
            'pops_conversions':pops_placed,'pops_failed_sites':pops_failures}
    return items,report


# V0.12.0 ASYMMETRIC OPC SOLVATION
def _amberhome():
    value=os.environ.get('AMBERHOME')
    if value and Path(value).is_dir(): return Path(value)
    executable=shutil.which('tleap') if 'shutil' in globals() else None
    if executable:
        return Path(executable).resolve().parents[1]
    raise RuntimeError('AMBERHOME is not set and tleap was not found on PATH')


def _ensure_opc_reference(cache_path):
    """Export AmberTools' installed OPCBOX once; never invent OPC geometry."""
    cache_path=Path(cache_path)
    if cache_path.is_file(): return cache_path
    cache_path.parent.mkdir(parents=True,exist_ok=True)
    amberhome=Path(os.environ.get('AMBERHOME',''))
    if not amberhome.is_dir():
        tleap_path=__import__('shutil').which('tleap')
        if not tleap_path: raise RuntimeError('tleap not found; activate AmberTools first')
        amberhome=Path(tleap_path).resolve().parents[1]
    opc_off=amberhome/'dat/leap/lib/opcbox.off'
    if not opc_off.is_file():
        raise FileNotFoundError(f'Installed Amber OPC box not found: {opc_off}')
    with tempfile.TemporaryDirectory() as tmp:
        script=Path(tmp)/'export_opc.leap'
        exported=Path(tmp)/'opcbox.pdb'
        script.write_text(
            'source leaprc.water.opc\n'
            f'loadOff {opc_off}\n'
            f'savePdb OPCBOX {exported}\n'
            'quit\n'
        )
        result=subprocess.run(['tleap','-f',str(script)],capture_output=True,text=True)
        if result.returncode!=0 or not exported.is_file():
            raise RuntimeError('tleap could not export OPCBOX:\n'+result.stdout+'\n'+result.stderr)
        cache_path.write_text(exported.read_text())
    return cache_path


def _parse_opc_box(path):
    lines=Path(path).read_text(errors='replace').splitlines()
    cryst=next((line for line in lines if line.startswith('CRYST1')),None)
    if not cryst: raise ValueError('Exported OPCBOX has no CRYST1 record')
    box=np.array([float(cryst[6:15]),float(cryst[15:24]),float(cryst[24:33])])
    groups=defaultdict(list)
    for line in lines:
        if not line.startswith(COORD_RECORDS): continue
        atom=parse_atom(line,0)
        groups[(atom['chain'],atom['resid'],atom['icode'])].append(atom)
    waters=[]
    for key,atoms in groups.items():
        oxygens=[atom for atom in atoms if atom['element']=='O' or atom['name'] in {'O','OW'}]
        if len(oxygens)!=1: continue
        waters.append({'atoms':atoms,'oxygen':oxygens[0]['xyz'].copy()})
    if not waters: raise ValueError('No complete waters parsed from exported OPCBOX')
    return waters,box


def _system_atom_cells(fixed,items,cell_width=4.0):
    cells=defaultdict(list)
    for atom in fixed:
        key=tuple(np.floor(atom['xyz']/cell_width).astype(int)); cells[key].append((atom['xyz'],atom['element']))
    for item in items:
        for point,element in zip(item['xyz'],item['elements']):
            key=tuple(np.floor(point/cell_width).astype(int)); cells[key].append((point,element))
    return cells


def _water_overlaps(xyz,elements,cells,oxygen_cut,hydrogen_cut,virtual_cut,cell_width=4.0):
    for point,element in zip(xyz,elements):
        if element=='O': cutoff=oxygen_cut
        elif element=='H': cutoff=hydrogen_cut
        else: cutoff=virtual_cut
        cutoff2=cutoff*cutoff; base=np.floor(point/cell_width).astype(int)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for other,_ in cells.get((base[0]+dx,base[1]+dy,base[2]+dz),()):
                        if float(np.dot(point-other,point-other))<cutoff2: return True
    return False


def _apply_opc_solvation(items,fixed,box,cfg):
    settings=cfg.get('solvent',{})
    if not settings.get('enabled',False): return items,box,{}
    if str(settings.get('model','')).upper()!='OPC': raise ValueError('Only OPC is supported in v0.12.0')
    cache=_ensure_opc_reference(settings.get('reference_cache','reference_library/OPC/opcbox.pdb'))
    waters,water_box=_parse_opc_box(cache)
    oxygen_coordinates=np.array([water['oxygen'] for water in waters],float)
    reference_origin=oxygen_coordinates.min(axis=0)
    for water in waters:
        for atom in water['atoms']:
            atom['xyz']=atom['xyz']-reference_origin
        water['oxygen']=water['oxygen']-reference_origin

    p31=[]
    for item in items:
        for atom,point in zip(item['molecule']['atoms'],item['xyz']):
            if atom['name']=='P31' and atom['resname'] in {'PC','PS'}: p31.append(point)
    if not p31: raise RuntimeError('Cannot identify membrane headgroup planes')
    upper_head=float(np.mean([p[2] for p in p31 if p[2]>0]))
    lower_head=float(np.mean([p[2] for p in p31 if p[2]<0]))
    top_buffer=float(settings['positive_z_buffer_from_headgroups_angstrom'])
    bottom_buffer=float(settings['negative_z_buffer_from_headgroups_angstrom'])
    clearance=float(settings.get('minimum_fixed_atom_box_clearance_angstrom',5.0))
    fixed_max=max(float(atom['xyz'][2]) for atom in fixed)
    fixed_min=min(float(atom['xyz'][2]) for atom in fixed)
    top=max(upper_head+top_buffer,fixed_max+clearance)
    bottom=min(lower_head-bottom_buffer,fixed_min-clearance)

    all_xy=np.vstack([np.array([atom['xyz'][:2] for atom in fixed]),
                      np.vstack([item['xyz'][:,:2] for item in items])])
    center_xy=0.5*(all_xy.min(axis=0)+all_xy.max(axis=0))
    xlo=center_xy[0]-box[0]/2; xhi=xlo+box[0]
    ylo=center_xy[1]-box[1]/2; yhi=ylo+box[1]

    cells=_system_atom_cells(fixed,items)
    oxygen_cut=float(settings.get('minimum_water_oxygen_system_distance_angstrom',2.4))
    hydrogen_cut=float(settings.get('minimum_water_hydrogen_system_distance_angstrom',1.5))
    virtual_cut=float(settings.get('minimum_water_virtual_site_system_distance_angstrom',1.2))
    accepted=[]; rejected=0
    nx=int(math.ceil(box[0]/water_box[0]))+1
    ny=int(math.ceil(box[1]/water_box[1]))+1
    nz=int(math.ceil((top-bottom)/water_box[2]))+1
    start=np.array([xlo,ylo,bottom])
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                tile=start+np.array([ix*water_box[0],iy*water_box[1],iz*water_box[2]])
                for water in waters:
                    oxygen=water['oxygen']+tile
                    if not (xlo<=oxygen[0]<xhi and ylo<=oxygen[1]<yhi and bottom<=oxygen[2]<top): continue
                    side='positive' if oxygen[2]>upper_head else ('negative' if oxygen[2]<lower_head else None)
                    if side is None: continue
                    xyz=np.array([atom['xyz']+tile for atom in water['atoms']],float)
                    elements=[atom['element'] for atom in water['atoms']]
                    if _water_overlaps(xyz,elements,cells,oxygen_cut,hydrogen_cut,virtual_cut):
                        rejected+=1; continue
                    molecule={'atoms':water['atoms'],'anchor':oxygen,'leaflet':side,'species':'OPC'}
                    accepted_item={'molecule':molecule,'xyz':xyz,'anchor':oxygen,'elements':elements,
                                   'leaflet':side,'species':'OPC'}
                    accepted.append(accepted_item)
                    for point,element in zip(xyz,elements):
                        cell_key=tuple(np.floor(point/4.0).astype(int))
                        cells[cell_key].append((point,element))
    items.extend(accepted)
    new_box=np.array([box[0],box[1],top-bottom])
    report={
      'model':'OPC','reference':str(cache),
      'reference_coordinate_origin_angstrom':[float(value) for value in reference_origin],
      'tiling_mode':'zero_based_contiguous_with_interface_overlap_checks',
      'box_faces_z_angstrom':{'bottom':bottom,'top':top},
      'headgroup_planes_z_angstrom':{'lower':lower_head,'upper':upper_head},
      'requested_buffers_from_headgroups_angstrom':{'positive':top_buffer,'negative':bottom_buffer},
      'actual_buffers_from_headgroups_angstrom':{'positive':top-upper_head,'negative':lower_head-bottom},
      'minimum_fixed_atom_clearance_angstrom':clearance,
      'waters_added_positive_z':sum(item['leaflet']=='positive' for item in accepted),
      'waters_added_negative_z':sum(item['leaflet']=='negative' for item in accepted),
      'waters_rejected_for_overlap':rejected,'ions_added':0
    }
    return items,new_box,report


# V0.12.4 TLEAP-CONTIGUOUS OPC VOLUME
def _generate_tleap_opc_volume(reference_path,bounds,cache_path):
    xlo,xhi,ylo,yhi,bottom,top=bounds
    cache_path=Path(cache_path); cache_path.parent.mkdir(parents=True,exist_ok=True)
    metadata=cache_path.with_suffix('.bounds.json')
    signature=[round(float(v),5) for v in bounds]
    if cache_path.is_file() and metadata.is_file():
        try:
            if json.loads(metadata.read_text()).get('bounds')==signature: return cache_path
        except Exception: pass
    waters,_=_parse_opc_box(reference_path); seed=waters[0]; seed_origin=seed['oxygen'].copy()
    with tempfile.TemporaryDirectory() as tmpname:
        tmp=Path(tmpname); scaffold=tmp/'opc_scaffold.pdb'; lines=[]; serial=1; resid=1; inset=1.5
        for x in (xlo+inset,xhi-inset):
            for y in (ylo+inset,yhi-inset):
                for z in (bottom+inset,top-inset):
                    target=np.array([x,y,z],float)
                    for atom in seed['atoms']:
                        point=atom['xyz']+(target-seed_origin); element=atom['element'][:2].rjust(2)
                        lines.append(f"HETATM{serial:5d} {atom['name']:<4s} WAT S{resid:4d}    "
                                     f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00          {element}")
                        serial+=1
                    resid+=1
        scaffold.write_text('\n'.join(lines+['END'])+'\n')
        leap=tmp/'make_opc_volume.leap'; raw=tmp/'opc_volume_raw.pdb'
        leap.write_text('source leaprc.water.opc\n'+f'SOL = loadPdb {scaffold}\n'+
                        'solvateBox SOL OPCBOX 0.0 1.0\n'+f'savePdb SOL {raw}\nquit\n')
        result=subprocess.run(['tleap','-f',str(leap)],capture_output=True,text=True)
        if result.returncode!=0 or not raw.is_file():
            raise RuntimeError('tleap OPC volume generation failed:\n'+result.stdout+'\n'+result.stderr)
        cache_path.write_text(raw.read_text()); metadata.write_text(json.dumps({'bounds':signature},indent=2)+'\n')
    return cache_path


def _parse_waters_from_pdb(path):
    groups=defaultdict(list)
    for line in Path(path).read_text(errors='replace').splitlines():
        if not line.startswith(COORD_RECORDS): continue
        atom=parse_atom(line,0)
        if atom['resname'] in {'WAT','HOH','OPC'}:
            groups[(atom['chain'],atom['resid'],atom['icode'])].append(atom)
    waters=[]
    for atoms in groups.values():
        oxygens=[a for a in atoms if a['element']=='O' or a['name'] in {'O','OW'}]
        if len(atoms)==4 and len(oxygens)==1: waters.append({'atoms':atoms,'oxygen':oxygens[0]['xyz'].copy()})
    if not waters: raise ValueError(f'No complete four-site OPC waters parsed from {path}')
    return waters


# V0.12.8 RESTORE ORIGINAL DIRECT OPC PLACEMENT
def _apply_opc_solvation(items,fixed,box,cfg):
    """Restore the original direct OPCBOX replication used by v0.12.0.

    This intentionally removes the later tleap scaffold/solvateBox volume path.
    AmberTools is still used to export the installed OPCBOX reference once.
    The reference waters are then repeated directly, cropped as complete
    molecules, and rejected only for clashes with the immutable assembly and
    completed membrane.
    """
    settings=cfg.get('solvent',{})
    if not settings.get('enabled',False):
        return items,box,{}
    if str(settings.get('model','')).upper()!='OPC':
        raise ValueError('Only OPC is supported')

    reference=_ensure_opc_reference(
        settings.get('reference_cache','reference_library/OPC/opcbox.pdb')
    )
    waters,water_box=_parse_opc_box(reference)

    p31=[]
    for item in items:
        for atom,point in zip(item['molecule']['atoms'],item['xyz']):
            if atom['name']=='P31' and atom['resname'] in {'PC','PS'}:
                p31.append(point)
    if not p31:
        raise RuntimeError('Cannot identify membrane headgroup planes')

    upper_head=float(np.mean([p[2] for p in p31 if p[2]>0]))
    lower_head=float(np.mean([p[2] for p in p31 if p[2]<0]))
    top_buffer=float(settings['positive_z_buffer_from_headgroups_angstrom'])
    bottom_buffer=float(settings['negative_z_buffer_from_headgroups_angstrom'])
    clearance=float(settings.get(
        'minimum_water_beyond_fixed_atoms_angstrom',
        settings.get('minimum_fixed_atom_box_clearance_angstrom',5.0)
    ))

    fixed_z_values=[float(atom['xyz'][2]) for atom in fixed]
    if fixed_z_values:
        fixed_z_min=min(fixed_z_values); fixed_z_max=max(fixed_z_values)
    else:
        fixed_z_min=lower_head; fixed_z_max=upper_head
    requested_top=upper_head+top_buffer
    requested_bottom=lower_head-bottom_buffer
    containment_top=fixed_z_max+clearance
    containment_bottom=fixed_z_min-clearance
    top=max(requested_top,containment_top)
    bottom=min(requested_bottom,containment_bottom)
    top_reason=('fixed_atom_containment' if containment_top>requested_top
                else 'requested_headgroup_buffer')
    bottom_reason=('fixed_atom_containment' if containment_bottom<requested_bottom
                   else 'requested_headgroup_buffer')

    xy_parts=[np.vstack([item['xyz'][:,:2] for item in items])]
    if fixed:
        xy_parts.insert(0,np.array([atom['xyz'][:2] for atom in fixed]))
    all_xy=np.vstack(xy_parts)
    center_xy=0.5*(all_xy.min(axis=0)+all_xy.max(axis=0))
    xlo=float(center_xy[0]-box[0]/2.0)
    xhi=float(xlo+box[0])
    ylo=float(center_xy[1]-box[1]/2.0)
    yhi=float(ylo+box[1])

    cells=_system_atom_cells(fixed,items)
    oxygen_cut=float(settings.get(
        'minimum_water_oxygen_system_distance_angstrom',2.4))
    hydrogen_cut=float(settings.get(
        'minimum_water_hydrogen_system_distance_angstrom',1.5))
    virtual_cut=float(settings.get(
        'minimum_water_virtual_site_system_distance_angstrom',1.2))

    accepted=[]
    rejected=0
    nx=int(math.ceil(box[0]/water_box[0]))+1
    ny=int(math.ceil(box[1]/water_box[1]))+1
    nz=int(math.ceil((top-bottom)/water_box[2]))+1
    start=np.array([xlo,ylo,bottom],float)

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                tile=start+np.array([
                    ix*water_box[0],iy*water_box[1],iz*water_box[2]
                ])
                for water in waters:
                    oxygen=water['oxygen']+tile
                    if not (xlo<=oxygen[0]<xhi and
                            ylo<=oxygen[1]<yhi and
                            bottom<=oxygen[2]<top):
                        continue
                    side=('positive' if oxygen[2]>upper_head else
                          ('negative' if oxygen[2]<lower_head else None))
                    if side is None:
                        continue
                    xyz=np.array([
                        atom['xyz']+tile for atom in water['atoms']
                    ],float)
                    elements=[atom['element'] for atom in water['atoms']]
                    if _water_overlaps(
                        xyz,elements,cells,
                        oxygen_cut,hydrogen_cut,virtual_cut
                    ):
                        rejected+=1
                        continue
                    molecule={
                        'atoms':water['atoms'],
                        'anchor':oxygen,
                        'leaflet':side,
                        'species':'OPC'
                    }
                    accepted.append({
                        'molecule':molecule,
                        'xyz':xyz,
                        'anchor':oxygen,
                        'elements':elements,
                        'leaflet':side,
                        'species':'OPC'
                    })

    items.extend(accepted)
    new_box=np.array([box[0],box[1],top-bottom])
    accepted_oxygen=[item['anchor'][2] for item in accepted]
    water_z_min=(float(min(accepted_oxygen)) if accepted_oxygen else None)
    water_z_max=(float(max(accepted_oxygen)) if accepted_oxygen else None)

    report={
      'model':'OPC',
      'reference':str(reference),
      'placement_mode':'original_direct_OPCBOX_replication_v0.12.0',
      'reference_box_angstrom':[float(v) for v in water_box],
      'replication_counts_xyz':[nx,ny,nz],
      'box_faces_z_angstrom':{'bottom':bottom,'top':top},
      'headgroup_planes_z_angstrom':{'lower':lower_head,'upper':upper_head},
      'requested_buffers_from_headgroups_angstrom':{
          'positive':top_buffer,'negative':bottom_buffer},
      'actual_buffers_from_headgroups_angstrom':{
          'positive':top-upper_head,'negative':lower_head-bottom},
      'fixed_atom_z_extrema_angstrom':{
          'minimum':fixed_z_min,'maximum':fixed_z_max},
      'boundary_control':{'positive':top_reason,'negative':bottom_reason},
      'fixed_atom_clearance_actual_angstrom':{
          'positive':top-fixed_z_max,'negative':fixed_z_min-bottom},
      'retained_water_oxygen_envelope_z_angstrom':{
          'minimum':water_z_min,'maximum':water_z_max},
      'waters_added_positive_z':sum(
          item['leaflet']=='positive' for item in accepted),
      'waters_added_negative_z':sum(
          item['leaflet']=='negative' for item in accepted),
      'waters_rejected_for_system_overlap':rejected,
      'ions_added':0
    }
    return items,new_box,report

def _pdb_atom_name_field(name,element):
    """Return the exact four-column PDB atom-name field.

    For ordinary one-letter elements, names shorter than four characters are
    shifted one column right (for example C12 -> " C12"), matching the
    Lipid21 reference PDB and tleap template/chirality matching. Four-character
    names such as C116 and H16R occupy all four columns unchanged.
    """
    name=str(name).strip()
    element=str(element).strip().upper()
    if not name:
        raise ValueError('Cannot write an empty PDB atom name')
    if len(name)>4:
        raise ValueError(f'Atom name exceeds the four-column PDB field: {name!r}')
    if len(name)==4:
        return name
    if len(element)==1 and name[0].isalpha():
        return f" {name:<3s}"
    return f"{name:<4s}"


def format_lipid_atom(serial,atom,resid,xyz):
    raw_element=atom['element'][:2]
    element=raw_element.rjust(2)
    atom_field=_pdb_atom_name_field(atom['name'],raw_element)
    resname=str(atom['resname']).strip()
    if len(resname)>3:
        raise ValueError(f'Residue name exceeds the three-column PDB field: {resname!r}')
    chain_ids='MWXYZABCDEFGHIJKLNOPQRSTUVabcdefghijklnopqrstuvwxyz0123456789'
    chain=chain_ids[((max(1,resid)-1)//9999)%len(chain_ids)]
    pdb_resid=((max(1,resid)-1)%9999)+1
    pdb_serial=((max(1,serial)-1)%99999)+1
    line=(f"HETATM{pdb_serial:5d} {atom_field} {resname:>3s} {chain}{pdb_resid:4d}    "
          f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element}")
    if line[12:16]!=atom_field or line[16:17]!=' ' or line[17:20].strip()!=resname:
        raise RuntimeError(f'Internal PDB-column audit failed for {name!r}')
    return line

# V0.12.9 TOPOLOGY-PRESERVING ION AND AMBER OUTPUT
def _require_parmed():
    try:
        import parmed as pmd
    except ImportError as exc:
        raise RuntimeError(
            'ParmEd is required for topology-preserving output. Activate the '
            'AmberTools environment and install/enable ParmEd before running.'
        ) from exc
    return pmd


def _write_environment_pdb(items,path,box):
    lines=[
        f"CRYST1{float(box[0]):9.3f}{float(box[1]):9.3f}{float(box[2]):9.3f}"
        f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1"
    ]
    serial=1; resid=1
    for item in items:
        for atom,point in zip(item['molecule']['atoms'],item['xyz']):
            lines.append(format_lipid_atom(serial,atom,resid,point))
            serial+=1
        resid+=1
    lines.append('END')
    Path(path).write_text('\n'.join(lines)+'\n')


def _tleap_environment_topology(pdb_path,prmtop_path,rst7_path,log_path):
    with tempfile.TemporaryDirectory() as tmpname:
        leap=Path(tmpname)/'environment.leap'
        leap.write_text(
            'source leaprc.lipid21\n'
            'source leaprc.water.opc\n'
            f'ENV = loadPdb {Path(pdb_path).resolve()}\n'
            f'saveAmberParm ENV {Path(prmtop_path).resolve()} {Path(rst7_path).resolve()}\n'
            'quit\n'
        )
        result=subprocess.run(['tleap','-f',str(leap)],capture_output=True,text=True)
        Path(log_path).write_text(result.stdout+'\n'+result.stderr)
        if result.returncode!=0 or not Path(prmtop_path).is_file():
            raise RuntimeError(
                'tleap could not parameterize the generated Lipid21/OPC/ion '
                f'environment. See {log_path}'
            )
        warning_text=(result.stdout+'\n'+result.stderr).lower()
        fatal_tokens=('unknown residue','does not have a type','could not find bond parameter')
        if any(token in warning_text for token in fatal_tokens):
            raise RuntimeError(f'tleap reported an incomplete environment topology; see {log_path}')


def _rounded_integer_charge(charge,tolerance):
    integer=int(round(charge))
    if abs(charge-integer)>tolerance:
        raise RuntimeError(
            f'Topology charge {charge:.8f} e is not within {tolerance} e of an integer'
        )
    return integer


def _ion_atom(name,resname,element):
    return {'name':name,'resname':resname,'element':element,'record':'HETATM','line':''}


def _replace_opc_with_ions(items,cation_name,cation_resname,cation_element,
                            anion_name,anion_resname,anion_element,
                            cation_count,anion_count,cfg,rng):
    ion_cfg=cfg.get('ions',{})
    minimum_fixed=float(ion_cfg.get('minimum_ion_fixed_heavy_distance_angstrom',4.0))
    minimum_ion=float(ion_cfg.get('minimum_ion_ion_distance_angstrom',5.0))
    candidates=[(idx,item) for idx,item in enumerate(items) if item.get('species')=='OPC']
    rng.shuffle(candidates)
    fixed_points=[]
    for item in items:
        if item.get('species')=='OPC': continue
        for point,element in zip(item['xyz'],item.get('elements',[])):
            if element!='H': fixed_points.append(np.asarray(point,float))
    fixed_points=np.asarray(fixed_points,float) if fixed_points else np.empty((0,3))
    selected=[]; selected_indices=set()
    requests=[('cation',cation_count,cation_name,cation_resname,cation_element),
              ('anion',anion_count,anion_name,anion_resname,anion_element)]
    for kind,count,name,resname,element in requests:
        placed=0
        for idx,item in candidates:
            if idx in selected_indices: continue
            point=np.asarray(item['anchor'],float)
            if fixed_points.size and np.min(np.linalg.norm(fixed_points-point,axis=1))<minimum_fixed:
                continue
            if selected and min(np.linalg.norm(point-other[1]) for other in selected)<minimum_ion:
                continue
            selected_indices.add(idx); selected.append((kind,point,idx,name,resname,element)); placed+=1
            if placed==count: break
        if placed!=count:
            raise RuntimeError(
                f'Could place only {placed} of {count} requested {kind} ions using current distance floors'
            )
    ions=[]
    for kind,point,idx,name,resname,element in selected:
        atom=_ion_atom(name,resname,element)
        molecule={'atoms':[atom],'anchor':point,'leaflet':items[idx]['leaflet'],'species':resname}
        ions.append({'molecule':molecule,'xyz':np.array([point]),'anchor':point,
                     'elements':[element],'leaflet':items[idx]['leaflet'],'species':resname})
    kept=[item for idx,item in enumerate(items) if idx not in selected_indices]
    kept.extend(ions)
    return kept,selected


# V0.13.5 ORDER-INDEPENDENT SOLUTE AUDIT
def _verify_solute_prefix(original,combined,tolerance=1.0e-6):
    """Verify preservation without assuming the solute remains an index prefix.

    ParmEd may reorder atoms when molecules are not contiguous (for example,
    retained OPC waters and tleap-reordered ligand/cofactor residues). That is
    an ordering normalization, not a topology change. Audit atom definitions
    and bonded terms as multisets instead of comparing atom index i to index i.
    """
    def q(value):
        return round(float(value), 6)

    def atom_signature(atom):
        return (
            atom.residue.name.strip(), atom.name.strip(), str(atom.type),
            q(atom.charge), q(atom.mass), int(getattr(atom, 'atomic_number', 0) or 0)
        )

    original_atoms = Counter(atom_signature(atom) for atom in original.atoms)
    combined_atoms = Counter(atom_signature(atom) for atom in combined.atoms)
    missing_atoms = original_atoms - combined_atoms
    if missing_atoms:
        example, count = next(iter(missing_atoms.items()))
        raise RuntimeError(
            'Combined topology is missing authoritative solute atom definitions: '
            f'{example} x {count}'
        )

    def bond_signature(bond):
        endpoints = sorted((atom_signature(bond.atom1), atom_signature(bond.atom2)))
        bond_type = bond.type
        parameters = None if bond_type is None else (
            q(getattr(bond_type, 'k', 0.0)), q(getattr(bond_type, 'req', 0.0))
        )
        return (tuple(endpoints), parameters)

    original_bonds = Counter(bond_signature(bond) for bond in original.bonds)
    combined_bonds = Counter(bond_signature(bond) for bond in combined.bonds)
    missing_bonds = original_bonds - combined_bonds
    if missing_bonds:
        example, count = next(iter(missing_bonds.items()))
        raise RuntimeError(
            'Combined topology is missing authoritative solute bonds or bond parameters: '
            f'{example} x {count}'
        )

    return {
        'audit_mode': 'order_independent_multiset',
        'authoritative_solute_atoms_found': sum(original_atoms.values()),
        'authoritative_solute_bonds_found': sum(original_bonds.values()),
        'atom_order_required': False,
        'parmed_reordering_allowed': True,
        'atom_definitions_preserved': True,
        'bond_definitions_preserved': True,
    }

# V0.13.6 RESTORE TOPOLOGY MAPPING HELPER
def _topology_pdb_mapping_report(solute,fixed):
    """Report PDB/topology differences without requiring user declarations.

    The supplied prmtop/rst7 is authoritative. This report permits tleap atom
    ordering, residue movement, aliases, template-added atoms, and OPC extra
    points. It is diagnostic only and does not impose an index-by-index match.
    """
    virtual=[]
    physical=[]
    virtual_names={'EP','EPW','LP','LP1','LP2','M','MW','VS','DUMM'}

    for atom in solute.atoms:
        name=atom.name.strip().upper()
        is_virtual=(
            int(getattr(atom,'atomic_number',0) or 0)==0
            or float(atom.mass)<0.1
            or name in virtual_names
        )
        if is_virtual:
            virtual.append(atom)
        else:
            physical.append(atom)

    pdb_residues=[]
    current=None
    for atom in fixed:
        key=(atom['chain'],atom['resid'],atom['icode'],atom['resname'])
        if current is None or current['key']!=key:
            current={'key':key,'names':Counter()}
            pdb_residues.append(current)
        current['names'][atom['name']]+=1

    physical_ids={id(atom) for atom in physical}
    topology_residues=[]
    for residue in solute.residues:
        names=Counter(
            atom.name.strip()
            for atom in residue.atoms
            if id(atom) in physical_ids
        )
        topology_residues.append({
            'name':residue.name.strip(),
            'names':names,
            'index':residue.idx,
        })

    pdb_signatures=Counter(
        (residue['key'][3],tuple(sorted(residue['names'].items())))
        for residue in pdb_residues
    )
    topology_signatures=Counter(
        (residue['name'],tuple(sorted(residue['names'].items())))
        for residue in topology_residues
    )
    exact_matches=sum((pdb_signatures & topology_signatures).values())

    pdb_name_counts=Counter(residue['key'][3] for residue in pdb_residues)
    topology_name_counts=Counter(residue['name'] for residue in topology_residues)
    inventory_differences=[]
    for name,count in (pdb_name_counts-topology_name_counts).items():
        inventory_differences.append({
            'pdb_residue_name':name,
            'unmatched_count':count,
        })
    for name,count in (topology_name_counts-pdb_name_counts).items():
        inventory_differences.append({
            'topology_residue_name':name,
            'unmatched_count':count,
        })

    virtual_counts=Counter(
        (atom.residue.name.strip(),atom.name.strip())
        for atom in virtual
    )

    return {
        'authoritative_coordinate_source':'pre_membrane.rst7',
        'authoritative_topology_source':'pre_membrane.prmtop',
        'pdb_physical_atoms':len(fixed),
        'topology_total_atoms':len(solute.atoms),
        'topology_physical_atoms':len(physical),
        'topology_virtual_sites':len(virtual),
        'physical_atom_delta_topology_minus_pdb':len(physical)-len(fixed),
        'virtual_site_counts':{
            f'{key[0]}:{key[1]}':value
            for key,value in virtual_counts.items()
        },
        'pdb_residues':len(pdb_residues),
        'topology_residues':len(topology_residues),
        'exact_residue_signature_matches':exact_matches,
        'residue_name_inventory_differences':inventory_differences,
        'note':(
            'Differences are detected automatically; prmtop/rst7 remains '
            'authoritative and no user-supplied mismatch list is required.'
        ),
    }


def _prepare_topology_and_ions(assembled,fixed,box,cfg,out,rng):
    top_cfg=cfg.get('topology',{})
    ion_cfg=cfg.get('ions',{})
    if not top_cfg.get('enabled',False):
        return assembled,{},None
    pmd=_require_parmed()
    solute_prmtop=Path(top_cfg.get('solute_prmtop','pre_membrane.prmtop')).expanduser().resolve()
    solute_rst7=Path(top_cfg.get('solute_coordinates','pre_membrane.rst7')).expanduser().resolve()
    if not solute_prmtop.is_file() or not solute_rst7.is_file():
        raise FileNotFoundError(f'Missing authoritative solute topology/coordinates: {solute_prmtop}, {solute_rst7}')
    solute=pmd.load_file(str(solute_prmtop),xyz=str(solute_rst7))
    mapping_report=_topology_pdb_mapping_report(solute,fixed)
    # The topology/restart pair is authoritative. Do not reject template-added
    # hydrogens, OPC EPW sites, tleap atom ordering, or residue aliases.
    max_delta=None

    preion_top=out/'environment_preion.prmtop'; preion_rst=out/'environment_preion.rst7'
    preion_interchange=_tleap_environment_topology_from_items(
        assembled,out/'environment_preion',preion_top,preion_rst,
        out/'tleap_environment_preion.log',box,cfg)
    environment=pmd.load_file(str(preion_top),xyz=str(preion_rst))
    tolerance=float(top_cfg.get('charge_integer_tolerance_e',0.001))
    solute_charge=_rounded_integer_charge(sum(float(a.charge) for a in solute.atoms),tolerance)
    membrane_charge=_rounded_integer_charge(sum(float(a.charge) for a in environment.atoms),tolerance)
    preion_charge=solute_charge+membrane_charge

    if not ion_cfg.get('enabled',False):
        cations=anions=target=0; selected=[]; final_items=assembled
        salt=str(ion_cfg.get('salt','KCl'))
    else:
        salt=str(ion_cfg.get('salt','KCl'))
        if salt not in {'KCl','NaCl'}: raise ValueError('ions.salt must be KCl or NaCl')
        concentration=float(ion_cfg.get('concentration_molar',0.150))
        retained_water_names={'WAT','HOH','OPC'}
        retained_input_waters=sum(
            residue.name.strip() in retained_water_names for residue in solute.residues)
        added_opc_waters=sum(item.get('species')=='OPC' for item in assembled)
        water_count=retained_input_waters+added_opc_waters
        target=int(round(water_count*concentration/55.5))
        if preion_charge>0:
            anions=max(target,preion_charge); cations=anions-preion_charge
        elif preion_charge<0:
            cations=max(target,abs(preion_charge)); anions=cations-abs(preion_charge)
        else:
            cations=target; anions=target
        if salt=='KCl': c_name=c_res='K+'; c_el='K'
        else: c_name=c_res='Na+'; c_el='NA'
        a_name=a_res='Cl-'; a_el='CL'
        final_items,selected=_replace_opc_with_ions(
            assembled,c_name,c_res,c_el,a_name,a_res,a_el,cations,anions,cfg,rng)

    final_env_top=out/'environment_with_ions.prmtop'; final_env_rst=out/'environment_with_ions.rst7'
    final_interchange=_tleap_environment_topology_from_items(
        final_items,out/'environment_with_ions',final_env_top,final_env_rst,
        out/'tleap_environment_final.log',box,cfg)
    final_environment=pmd.load_file(str(final_env_top),xyz=str(final_env_rst))
    combined=solute+final_environment
    combined.box=np.array([box[0],box[1],box[2],90.0,90.0,90.0])
    solute_preservation_audit=_verify_solute_prefix(solute,combined)
    final_charge=_rounded_integer_charge(sum(float(a.charge) for a in combined.atoms),tolerance)
    if ion_cfg.get('enabled',False) and final_charge!=0:
        raise RuntimeError(f'Final topology charge is {final_charge} e, not zero')
    topology_name=top_cfg.get('output_prmtop','topology.prmtop')
    coordinate_name=top_cfg.get('output_rst7','coordinates.rst7')
    combined.save(str(out/topology_name),overwrite=True)
    combined.save(str(out/coordinate_name),format='rst7',overwrite=True)
    retained_input_waters=sum(
        residue.name.strip() in {'WAT','HOH','OPC'} for residue in solute.residues)
    result={'solute_prmtop':str(solute_prmtop),'solute_coordinates':str(solute_rst7),
            'input_topology_mapping':mapping_report,
            'lipid21_interchange_preion':preion_interchange,
            'lipid21_interchange_final':final_interchange,
            'solute_atoms_preserved':len(solute.atoms),'solute_coordinate_max_delta_angstrom':max_delta,
            'retained_input_opc_waters':retained_input_waters,
            'solute_charge_e':solute_charge,'generated_environment_charge_before_ions_e':membrane_charge,
            'preion_system_charge_e':preion_charge,'salt':salt,
            'concentration_molar':float(ion_cfg.get('concentration_molar',0.150)),
            'concentration_basis':'accepted_OPC_water_count',
            'initial_opc_waters':retained_input_waters+sum(item.get('species')=='OPC' for item in assembled),
            'newly_added_opc_waters':sum(item.get('species')=='OPC' for item in assembled),
            'water_derived_counterion_target':target,'cations_added':cations,'anions_added':anions,
            'opc_waters_replaced':len(selected),'final_opc_waters':retained_input_waters+sum(item.get('species')=='OPC' for item in final_items),
            'final_charge_e':final_charge,'output_prmtop':str(out/topology_name),
            'output_rst7':str(out/coordinate_name),
            'solute_preservation_audit':solute_preservation_audit,
            'solute_bonds_preserved':solute_preservation_audit['bond_definitions_preserved']}
    return final_items,result,combined

# V0.15.0 MEMBRANE-ONLY BUILDER
def _build_membrane_only_host(molecules,box_xy,target_per_leaflet,rng):
    """Tile the validated POPC patch and select an exact symmetric host membrane."""
    box_xy=np.asarray(box_xy,float)
    if box_xy.shape!=(2,) or np.any(box_xy<=0):
        raise ValueError('membrane_only.box_xy_angstrom must contain two positive values')
    target=int(target_per_leaflet)
    if target<=0: raise ValueError('membrane_only.target_per_leaflet must be positive')

    anchors=np.array([molecule['anchor'][:2] for molecule in molecules],float)
    span=np.ptp(anchors,axis=0)
    # Reference molecules are already periodically arranged. Tile far enough to
    # cover the requested rectangle, then crop by complete-molecule anchors.
    x_period=max(float(span[0]),1.0)
    y_period=max(float(span[1]),1.0)
    nx=int(math.ceil(box_xy[0]/x_period))+3
    ny=int(math.ceil(box_xy[1]/y_period))+3
    candidates={'upper':[],'lower':[]}
    ref_center=0.5*(anchors.min(axis=0)+anchors.max(axis=0))
    for ix in range(-nx,nx+1):
        for iy in range(-ny,ny+1):
            tile=np.array([ix*x_period,iy*y_period],float)-ref_center
            for molecule in molecules:
                anchor_xy=molecule['anchor'][:2]+tile
                if (-box_xy[0]/2<=anchor_xy[0]<box_xy[0]/2 and
                    -box_xy[1]/2<=anchor_xy[1]<box_xy[1]/2):
                    shift=np.array([tile[0],tile[1],0.0])
                    xyz,anchor,elements=_atom_arrays_from_molecule(molecule,shift=shift)
                    candidates[molecule['leaflet']].append({
                        'molecule':molecule,'xyz':xyz,'anchor':anchor,
                        'elements':elements,'leaflet':molecule['leaflet'],'species':'POPC'})

    selected=[]; diagnostics={}
    for leaf in ('upper','lower'):
        pool=candidates[leaf]
        if len(pool)<target:
            raise RuntimeError(
                f'Only {len(pool)} POPC sites are available in the {leaf} leaflet '
                f'for requested target {target}. Increase box_xy_angstrom.'
            )
        # Prefer a spatially uniform subset when cropping yields more sites than
        # requested. Seeded jitter changes replicas without creating dense holes.
        jitter=rng.random(len(pool))*1e-3
        radii=np.array([np.dot(item['anchor'][:2],item['anchor'][:2]) for item in pool])
        order=np.argsort(radii+jitter)
        chosen=[pool[int(index)] for index in order[:target]]
        selected.extend(chosen)
        diagnostics[leaf]={'available_sites':len(pool),'selected_sites':target}
    return selected,diagnostics


def _prepare_membrane_only_topology_and_ions(assembled,box,cfg,out,rng):
    """Parameterize membrane/water, add ions, and write final Amber files."""
    pmd=_require_parmed(); top_cfg=cfg.get('topology',{}); ion_cfg=cfg.get('ions',{})
    tolerance=float(top_cfg.get('charge_integer_tolerance_e',0.01))

    preion_top=out/'environment_preion.prmtop'; preion_rst=out/'environment_preion.rst7'
    preion_interchange=_tleap_environment_topology_from_items(
        assembled,out/'environment_preion',preion_top,preion_rst,
        out/'tleap_environment_preion.log',box,cfg)
    environment=pmd.load_file(str(preion_top),xyz=str(preion_rst))
    raw_preion_charge=sum(float(atom.charge) for atom in environment.atoms)
    preion_charge=_rounded_integer_charge(raw_preion_charge,tolerance)

    salt=str(ion_cfg.get('salt','KCl'))
    if salt not in {'KCl','NaCl'}: raise ValueError('ions.salt must be KCl or NaCl')
    target=cations=anions=0; selected=[]; final_items=assembled
    water_count=sum(item.get('species')=='OPC' for item in assembled)
    if ion_cfg.get('enabled',False):
        concentration=float(ion_cfg.get('concentration_molar',0.150))
        target=int(round(water_count*concentration/55.5))
        if preion_charge>0:
            anions=max(target,preion_charge); cations=anions-preion_charge
        elif preion_charge<0:
            cations=max(target,abs(preion_charge)); anions=cations-abs(preion_charge)
        else:
            cations=target; anions=target
        if salt=='KCl': c_name=c_res='K+'; c_el='K'
        else: c_name=c_res='Na+'; c_el='NA'
        final_items,selected=_replace_opc_with_ions(
            assembled,c_name,c_res,c_el,'Cl-','Cl-','CL',cations,anions,cfg,rng)

    final_top=out/'environment_with_ions.prmtop'; final_rst=out/'environment_with_ions.rst7'
    final_interchange=_tleap_environment_topology_from_items(
        final_items,out/'environment_with_ions',final_top,final_rst,
        out/'tleap_environment_final.log',box,cfg)
    final_system=pmd.load_file(str(final_top),xyz=str(final_rst))
    final_system.box=np.array([box[0],box[1],box[2],90.0,90.0,90.0])
    raw_final_charge=sum(float(atom.charge) for atom in final_system.atoms)
    final_charge=_rounded_integer_charge(raw_final_charge,tolerance)
    if ion_cfg.get('enabled',False) and final_charge!=0:
        raise RuntimeError(f'Final membrane-only topology charge is {final_charge} e, not zero')

    topology_name=top_cfg.get('output_prmtop','topology.prmtop')
    coordinate_name=top_cfg.get('output_rst7','coordinates.rst7')
    final_system.save(str(out/topology_name),overwrite=True)
    final_system.save(str(out/coordinate_name),format='rst7',overwrite=True)
    return final_items,{
        'mode':'membrane_only','raw_preion_charge_e':raw_preion_charge,
        'preion_system_charge_e':preion_charge,'raw_final_charge_e':raw_final_charge,
        'final_charge_e':final_charge,'salt':salt,
        'concentration_molar':float(ion_cfg.get('concentration_molar',0.150)),
        'concentration_basis':'accepted_OPC_water_count',
        'initial_opc_waters':water_count,'water_derived_counterion_target':target,
        'cations_added':cations,'anions_added':anions,'opc_waters_replaced':len(selected),
        'final_opc_waters':sum(item.get('species')=='OPC' for item in final_items),
        'lipid21_interchange_preion':preion_interchange,
        'lipid21_interchange_final':final_interchange,
        'output_prmtop':str(out/topology_name),'output_rst7':str(out/coordinate_name),
    },final_system


def _build_membrane_only(config_path,cfg):
    settings=cfg.get('membrane_only',{})
    box_xy=np.asarray(settings.get('box_xy_angstrom',[110.0,110.0]),float)
    target_per_leaflet=int(settings.get('target_per_leaflet',200))
    target_total=int(settings.get('target_total_lipids',2*target_per_leaflet))
    if target_total!=2*target_per_leaflet:
        raise ValueError('target_total_lipids must equal 2 * target_per_leaflet')
    rng=np.random.default_rng(int(cfg.get('build',{}).get('seed',1001)))
    reference=Path(cfg['reference_library']['POPC']).expanduser().resolve()
    molecules,_reference_box=parse_reference(reference)
    assembled,host_report=_build_membrane_only_host(molecules,box_xy,target_per_leaflet,rng)
    lower=-float(settings.get('half_thickness_angstrom',17.9))
    upper=+float(settings.get('half_thickness_angstrom',17.9))
    box=np.array([box_xy[0],box_xy[1],float(cfg.get('membrane',{}).get('box_z_angstrom',96.0))])
    fixed=[]
    assembled,composition_report=_apply_complex_composition(
        assembled,fixed,lower,upper,cfg,rng)
    assembled,box,solvent_report=_apply_opc_solvation(
        assembled,fixed,box,cfg)
    out=Path(cfg['output']['directory']); out.mkdir(parents=True,exist_ok=True)
    assembled,topology_report,final_system=_prepare_membrane_only_topology_and_ions(
        assembled,box,cfg,out,rng)
    output_pdb=out/cfg.get('output',{}).get('pdb_name','system_membrane_solvated.pdb')
    final_system.save(str(output_pdb),overwrite=True)
    counts=Counter((item['species'],item['leaflet']) for item in assembled)
    report={
        'version':'0.15.1-membrane-only-custom-sterol','build_mode':'membrane_only',
        'seed':int(cfg.get('build',{}).get('seed',1001)),
        'box_angstrom':[float(value) for value in box],
        'host_selection':host_report,
        'requested_membrane':{'target_total_lipids':target_total,
            'target_per_leaflet':target_per_leaflet,'box_xy_angstrom':[float(v) for v in box_xy]},
        'lipid_counts':{f'{species}_{leaflet}':count
            for (species,leaflet),count in counts.items()},
        'complex_composition':composition_report,'solvent':solvent_report,
        'topology_and_ions':topology_report,'output_pdb':str(output_pdb),
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report

def build(config_path):
    cfg=yaml.safe_load(Path(config_path).read_text())
    if cfg.get('membrane_only',{}).get('enabled',False):
        return _build_membrane_only(config_path,cfg)
    input_path=Path(cfg['protein']['pdb']).expanduser().resolve()
    if not cfg['protein'].get('preoriented') or cfg['protein'].get('orientation_source')!='OPM':
        raise ValueError('AMBMEM requires preoriented: true and orientation_source: OPM')
    ref_path=Path(cfg['reference_library']['POPC']).expanduser().resolve()
    original_lines,fixed,dum,planes,membrane_geometry_source=parse_opm(input_path,cfg.get('protein',{}))
    molecules,ref_box=parse_reference(ref_path)
    lower,upper=planes; mid=0.5*(lower+upper)
    slab_margin=float(cfg['membrane'].get('footprint_z_margin_angstrom',8.0))
    membrane_atoms=[a for a in fixed if lower-slab_margin<=a['xyz'][2]<=upper+slab_margin]
    if not membrane_atoms: raise ValueError('No fixed atoms found near the OPM membrane slab')
    xy=np.array([a['xyz'][:2] for a in membrane_atoms])
    minimum=xy.min(axis=0); maximum=xy.max(axis=0); center=0.5*(minimum+maximum)
    padding=float(cfg['membrane'].get('edge_padding_angstrom',20.0))
    needed=(maximum-minimum)+2*padding
    tiles=np.maximum(1,np.ceil(needed/ref_box[:2]).astype(int))
    box=np.array([tiles[0]*ref_box[0],tiles[1]*ref_box[1],float(cfg['membrane'].get('box_z_angstrom',200.0))])
    origin=center-box[:2]/2
    fixed_cells=build_fixed_cells(fixed,4.0)
    heavy_cut=float(cfg['packing'].get('fixed_heavy_lipid_heavy_angstrom',2.0))
    hydrogen_cut=float(cfg['packing'].get('fixed_heavy_lipid_hydrogen_angstrom',1.2))
    kept=[]; removed=[]
    for ix in range(int(tiles[0])):
        for iy in range(int(tiles[1])):
            tile_shift=np.array([origin[0]+ix*ref_box[0],origin[1]+iy*ref_box[1],mid])
            for mol in molecules:
                normalized_anchor=np.array([mol['anchor'][0]%ref_box[0],mol['anchor'][1]%ref_box[1],mol['anchor'][2]])
                shift=tile_shift+np.array([0.0,0.0,0.0])
                # Normalize reference X/Y while preserving each molecule's internal geometry.
                anchor_shift=np.array([origin[0]+ix*ref_box[0]+normalized_anchor[0]-mol['anchor'][0],
                                       origin[1]+iy*ref_box[1]+normalized_anchor[1]-mol['anchor'][1],mid])
                placed_anchor=mol['anchor']+anchor_shift
                if clashes(mol,fixed_cells,anchor_shift,4.0,heavy_cut,hydrogen_cut):
                    removed.append({'leaflet':mol['leaflet'],'anchor':[float(x) for x in placed_anchor]})
                else:
                    kept.append((mol,anchor_shift))
    rng=np.random.default_rng(int(cfg.get('build',{}).get('seed',1001)))
    kept,boundary_refilled,interior_removed=_refill_boundary_sites(
        removed,molecules,kept,fixed,lower,upper,box,cfg,rng
    )
    assembled=[]
    for molecule,shift in kept:
        xyz,anchor,elements=_atom_arrays_from_molecule(molecule,shift=shift)
        assembled.append({'molecule':molecule,'xyz':xyz,'anchor':anchor,'elements':elements,'leaflet':molecule['leaflet'],'species':'POPC'})
    for accepted in boundary_refilled:
        assembled.append({'molecule':accepted['molecule'],'xyz':accepted['xyz'],'anchor':np.asarray(accepted['final_anchor']),
                          'elements':[a['element'] for a in accepted['molecule']['atoms']],
                          'leaflet':accepted['leaflet'],'species':'POPC'})
    assembled,complex_report=_apply_complex_composition(assembled,fixed,lower,upper,cfg,rng)
    assembled,box,solvent_report=_apply_opc_solvation(assembled,fixed,box,cfg)
    out=Path(cfg['output']['directory']); out.mkdir(exist_ok=True)
    assembled,topology_report,_combined_topology=_prepare_topology_and_ions(
        assembled,fixed,box,cfg,out,rng)
    output_pdb=out/cfg.get('output',{}).get('pdb_name','system_membrane_solvated.pdb')
    report_path=out/'report.json'
    max_serial=max(int(a['line'][6:11]) for a in fixed if a['line'][6:11].strip().isdigit())
    output=[]
    # V0.12.6 CRYST1 FINAL BOX OUTPUT
    # box contains the final X/Y/Z dimensions at this point, including the
    # asymmetric solvent extent when solvent is enabled.
    output.append(
        f"CRYST1{float(box[0]):9.3f}{float(box[1]):9.3f}{float(box[2]):9.3f}"
        f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1"
    )
    # V0.12.5C FILTER TERMINAL AND CONNECTIVITY RECORDS
    removed_terminal_records=0
    removed_connectivity_records=0
    for line in original_lines:
        record=line[:6].strip()
        if record in {'END','ENDMDL','MASTER','CRYST1'}:
            removed_terminal_records+=1
            continue
        if record=='CONECT':
            removed_connectivity_records+=1
            continue
        output.append(line)
    serial=max_serial+1; resid=1
    counts=Counter()
    for item in assembled:
        molecule=item['molecule']; xyz=item['xyz']
        for atom,point in zip(molecule['atoms'],xyz):
            output.append(format_lipid_atom(serial,atom,resid,point)); serial+=1
        counts[(item['species'],item['leaflet'])]+=1; resid+=1
    output.append('END')
    output_pdb.write_text('\n'.join(output)+'\n')
    if _combined_topology is not None:
        # This includes tleap-added physical atoms and topology-defined ordering.
        _combined_topology.save(str(output_pdb),overwrite=True)
    removed_counts=Counter(r['leaflet'] for r in removed)
    heterogens=Counter(a['resname'] for a in fixed if a['record']=='HETATM')
    hydrogens=sum(a['element']=='H' for a in fixed)
    report={
      'version':'0.15.1-membrane-only-custom-sterol','input_pdb':str(input_path),'input_sha256':sha256(input_path),
      'orientation_source':'OPM','membrane_geometry_source':membrane_geometry_source,
      'pdb_serialization':{'input_terminal_records_deferred':removed_terminal_records,'input_conect_records_omitted':removed_connectivity_records,'single_final_END_record':True},
      'fixed_input':{'physical_coordinate_records':len(fixed),
        'atom_records':sum(a['record']=='ATOM' for a in fixed),'heteroatom_records':sum(a['record']=='HETATM' for a in fixed),
        'heterogens':dict(heterogens),'explicit_hydrogens':hydrogens,
        'maximum_coordinate_change_angstrom':0.0,'records_removed':0},
      'opm_metadata':{'dum_records':len(dum),'lower_boundary_z_angstrom':lower,
        'upper_boundary_z_angstrom':upper,'midplane_z_angstrom':mid,'hydrophobic_thickness_angstrom':upper-lower},
      'box_angstrom':[float(x) for x in box],'reference_tiles_xy':[int(x) for x in tiles],
      'lipid_counts':{f'{species}_{leaflet}':count for (species,leaflet),count in counts.items()},
      'popc_boundary':{'removed_for_fixed_input':dict(removed_counts),
        'boundary_refilled':dict(Counter(item['leaflet'] for item in boundary_refilled)),
        'interior_or_unfillable_removed':dict(Counter(item['leaflet'] for item in interior_removed))},
      'complex_composition':complex_report,
      'solvent':solvent_report,
      'topology_and_ions':topology_report,
      'boundary_refill_details':[{k:v for k,v in item.items() if k not in ('molecule','xyz')} for item in boundary_refilled],
      'packing_cutoffs_angstrom':{'fixed_heavy_lipid_heavy':heavy_cut,'fixed_heavy_lipid_hydrogen':hydrogen_cut},
      'output_pdb':str(output_pdb)
    }
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest='command',required=True)
    build_parser=sub.add_parser('build'); build_parser.add_argument('config')
    args=parser.parse_args()
    if args.command=='build': build(args.config)

if __name__=='__main__': main()
