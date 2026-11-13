from pathlib import Path
# V0.1.7.1 FRESH FINAL SOLUTE RELOAD
# V0.1.7 STRUCTURAL WATER CLASH POLICY

# V0.1.7 STRUCTURAL WATER CLASH POLICY

# Fixed internal safety thresholds. They are intentionally not YAML options.
_STRUCTURAL_WATER_OXYGEN_HEAVY_CUTOFF_ANGSTROM = 1.0
_STRUCTURAL_WATER_ANY_SITE_CUTOFF_ANGSTROM = 0.5
_STRUCTURAL_WATER_WARNING_CUTOFF_ANGSTROM = 1.5


def _structural_water_clash_audit(solute, cfg, out=None, apply_removal=True):
    """Run one unified audit and optionally remove the union once."""
    pmd = _require_parmed()
    from scipy.spatial import cKDTree

    policy = str(cfg.get('structural_waters', {}).get(
        'hard_clash_policy', 'remove')).strip().lower()
    if policy not in {'error', 'remove', 'keep'}:
        raise ValueError('structural_waters.hard_clash_policy must be error, remove, or keep')

    working = solute.copy(pmd.Structure)
    water_names = set(_STRUCTURAL_WATER_RESNAMES)
    xyz = np.asarray(working.coordinates, dtype=float)
    nonwater_all = [a.idx for a in working.atoms
                    if a.residue.name.strip() not in water_names]
    nonwater_heavy = [a.idx for a in working.atoms
                      if a.residue.name.strip() not in water_names
                      and int(getattr(a, 'atomic_number', 0) or 0) > 1
                      and float(a.mass) > 0.0]
    all_tree = cKDTree(xyz[nonwater_all]) if nonwater_all else None
    heavy_tree = cKDTree(xyz[nonwater_heavy]) if nonwater_heavy else None

    records = []
    rejected = []
    warnings_ = []
    accepted = []
    for residue in working.residues:
        if residue.name.strip() not in water_names:
            continue
        oxygen = next((a for a in residue.atoms
                       if int(getattr(a, 'atomic_number', 0) or 0) == 8
                       or a.name.strip() in {'O', 'OW'}), None)
        if oxygen is None:
            raise RuntimeError(f'Input water residue {residue.idx+1} has no oxygen')

        oh_dist = float('inf'); oh_partner = None
        if heavy_tree is not None:
            d, local = heavy_tree.query(xyz[oxygen.idx], k=1)
            oh_dist = float(d)
            p = working.atoms[nonwater_heavy[int(local)]]
            oh_partner = f'{p.residue.name.strip()}_{p.residue.idx+1}@{p.name.strip()}'

        water_indices = [a.idx for a in residue.atoms]
        aa_dist = float('inf'); water_site = None; aa_partner = None
        if all_tree is not None and water_indices:
            ds, locals_ = all_tree.query(xyz[water_indices], k=1)
            wi = int(np.argmin(ds))
            aa_dist = float(ds[wi])
            water_site = residue.atoms[wi].name.strip()
            p = working.atoms[nonwater_all[int(locals_[wi])]]
            aa_partner = f'{p.residue.name.strip()}_{p.residue.idx+1}@{p.name.strip()}'

        criteria = []
        if oh_dist < _STRUCTURAL_WATER_OXYGEN_HEAVY_CUTOFF_ANGSTROM:
            criteria.append('oxygen_to_nonwater_heavy')
        if aa_dist < _STRUCTURAL_WATER_ANY_SITE_CUTOFF_ANGSTROM:
            criteria.append('any_water_site_to_nonwater_atom')

        rid = residue.idx + 1
        if criteria:
            rejected.append(rid)
            classification = ('reject_multiple_criteria' if len(criteria) > 1
                              else 'reject_' + criteria[0])
        elif oh_dist < _STRUCTURAL_WATER_WARNING_CUTOFF_ANGSTROM:
            warnings_.append(rid); classification = 'warning_retain'
        else:
            accepted.append(rid); classification = 'accept'

        records.append({
            'input_residue_index_1based': rid,
            'residue_name': residue.name.strip(),
            'classification': classification,
            'rejection_criteria': criteria,
            'oxygen_to_nonwater_heavy_distance_angstrom': oh_dist,
            'oxygen_nearest_nonwater_heavy_atom': oh_partner,
            'minimum_any_water_site_to_nonwater_atom_distance_angstrom': aa_dist,
            'closest_water_site': water_site,
            'nearest_nonwater_atom_for_any_site': aa_partner,
        })

    rejected = sorted(set(rejected))
    audit = {
        'policy': policy,
        'cutoffs_are_internal_not_yaml_configurable': True,
        'oxygen_to_nonwater_heavy_cutoff_angstrom': 1.0,
        'any_water_site_to_nonwater_atom_cutoff_angstrom': 0.5,
        'warning_cutoff_angstrom': 1.5,
        'input_structural_waters': len(records),
        'rejected_structural_waters': len(rejected),
        'hard_clashing_structural_waters': len(rejected),
        'removed_structural_waters': 0,
        'retained_structural_waters': len(records),
        'warning_retained_structural_waters': len(warnings_),
        'accepted_structural_waters': len(accepted),
        'rejected_residue_indices_1based': rejected,
        'hard_clash_residue_indices_1based': rejected,
        'warning_residue_indices_1based': warnings_,
        'waters': records,
        'ion_replacement_eligible': False,
        'counted_in_total_water_for_ion_target': True,
    }

    if rejected and policy == 'error':
        if out is not None:
            Path(out).mkdir(parents=True, exist_ok=True)
            (Path(out)/'structural_water_cleanup.json').write_text(json.dumps(audit, indent=2)+'\n')
        raise RuntimeError(f'Input contains {len(rejected)} rejected structural waters')
    if rejected and policy == 'remove' and apply_removal:
        working.strip(':' + ','.join(map(str, rejected)))
        audit['removed_structural_waters'] = len(rejected)
        audit['retained_structural_waters'] = len(records) - len(rejected)
        audit['simulation_ready_allowed'] = True
    elif rejected and policy == 'keep':
        audit['kept_hard_clashes'] = True
        audit['simulation_ready_allowed'] = False
    else:
        audit['simulation_ready_allowed'] = True
    if out is not None:
        Path(out).mkdir(parents=True, exist_ok=True)
        (Path(out)/'structural_water_cleanup.json').write_text(json.dumps(audit, indent=2)+'\n')
    return working, audit


def _run_cpptraj_text(commands, log_path, description):
    result = subprocess.run(['cpptraj'], input='\n'.join(commands)+'\n',
                            text=True, capture_output=True)
    Path(log_path).write_text(result.stdout+'\n'+result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f'{description} failed; see {log_path}')


def _prepare_authoritative_solute_pair(solute_prmtop, solute_rst7, cfg, out):
    """Order first, audit once, remove one union mask, reorder, and reload."""
    pmd = _require_parmed()
    out = Path(out)
    work = out/'authoritative_solute_preprocessing'
    work.mkdir(parents=True, exist_ok=True)
    ordered_top = work/'01_ordered_original.prmtop'
    ordered_rst = work/'01_ordered_original.rst7'
    clean_top = work/'02_unified_clean.prmtop'
    clean_rst = work/'02_unified_clean.rst7'
    final_top = work/'03_clean_ordered.prmtop'
    final_rst = work/'03_clean_ordered.rst7'

    _run_cpptraj_text([
        f'parm {solute_prmtop}', f'trajin {solute_rst7}',
        f'fixatomorder parmout {ordered_top}',
        f'trajout {ordered_rst} restart', 'run', 'quit'],
        work/'01_fixatomorder.log', 'Initial solute fixatomorder')
    ordered = pmd.load_file(str(ordered_top), xyz=str(ordered_rst))
    _, audit = _structural_water_clash_audit(ordered, cfg, apply_removal=False)
    rejected = audit['rejected_residue_indices_1based']
    policy = audit['policy']

    if rejected and policy == 'error':
        (out/'structural_water_cleanup.json').write_text(json.dumps(audit, indent=2)+'\n')
        raise RuntimeError(f'Input contains {len(rejected)} rejected structural waters')
    if rejected and policy == 'remove':
        mask = ':' + ','.join(map(str, rejected))
        _run_cpptraj_text([
            f'parm {ordered_top}', f'trajin {ordered_rst}', f'strip {mask}',
            f'trajout {clean_rst} restart', 'run', 'quit'],
            work/'02_strip_coordinates.log', 'Unified coordinate removal')
        _run_cpptraj_text([
            f'parm {ordered_top}', f'parmstrip {mask}', f'parmwrite out {clean_top}',
            'run', 'quit'], work/'02_strip_topology.log', 'Unified topology removal')
        audit['removed_structural_waters'] = len(rejected)
        audit['retained_structural_waters'] = audit['input_structural_waters']-len(rejected)
        audit['simulation_ready_allowed'] = True
    else:
        clean_top.write_bytes(ordered_top.read_bytes())
        clean_rst.write_bytes(ordered_rst.read_bytes())

    _run_cpptraj_text([
        f'parm {clean_top}', f'trajin {clean_rst}',
        f'fixatomorder parmout {final_top}', f'trajout {final_rst} restart',
        'run', 'quit'], work/'03_final_fixatomorder.log', 'Final solute fixatomorder')
    prepared = pmd.load_file(str(final_top), xyz=str(final_rst))
    _, remaining = _structural_water_clash_audit(
        prepared, {'structural_waters': {'hard_clash_policy': 'keep'}},
        apply_removal=False)
    if policy == 'remove' and remaining['rejected_structural_waters']:
        raise RuntimeError('Unified cleanup left rejected waters: ' +
                           str(remaining['rejected_residue_indices_1based']))

    audit['preprocessing_order'] = [
        'fixatomorder_original', 'single_unified_water_audit',
        'single_union_removal', 'fixatomorder_cleaned', 'exact_pair_reload']
    audit['prepared_solute_topology'] = str(final_top)
    audit['prepared_solute_coordinates'] = str(final_rst)
    audit['prepared_solute_atoms'] = len(prepared.atoms)
    audit['remaining_rejected_structural_waters'] = remaining['rejected_structural_waters']
    (out/'structural_water_cleanup.json').write_text(json.dumps(audit, indent=2)+'\n')
    return prepared, audit, final_top, final_rst

# V0.1.6 PROTECTED STRUCTURAL WATER BLOCKS
import json
# V0.1.5 GLOBAL EXACT CHOLESTEROL PLACEMENT
import copy
# V0.1.4 UNIFIED POPS POPE HEADGROUP REPLACEMENT
# V0.16.3 CANONICAL REPLICA ATOM ORDER
# V0.16.2 COMMENSURATE PERIODIC BUILDER
# V0.16.0 CPPTRAJ EXACT-FINAL-PAIR GATE
from .final_validation import validate_final_system
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

    POPC is written PA->PC->OL, POPE as PA->PE->OL, and POPS as
    PA->PS->OL. A TER record is placed only after a complete lipid.
    Internal coordinates are unchanged.
    """
    prefix=Path(prefix)
    modular=[]
    bulk=[]
    for item in items:
        residue_names={str(atom['resname']).strip() for atom in item['molecule']['atoms']}
        if item.get('species') in {'POPC','POPE','POPS'} or residue_names & {'PA','PC','PE','PS','OL'}:
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
            species=str(
                item.get('species', '')
            ).strip()

            headgroup_by_species={
                'POPC':'PC',
                'POPE':'PE',
                'POPS':'PS',
            }

            headgroup=headgroup_by_species.get(
                species
            )

            if headgroup is None:
                raise RuntimeError(
                    f'No Lipid21 headgroup definition '
                    f'exists for species {species!r}'
                )

            module_rank={
                'PA':0,
                headgroup:1,
                'OL':2,
            }
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
        'module_order':{
            'POPC':['PA','PC','OL'],
            'POPE':['PA','PE','OL'],
            'POPS':['PA','PS','OL'],
        },
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
        candidate['canonical_site_id']=items[index].get('canonical_site_id',(index,))
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

# V0.1.6 PROTECTED STRUCTURAL WATER BLOCKS
_STRUCTURAL_WATER_RESNAMES = ('WAT','HOH','OPC','TIP3','TIP4','TIP4P')


def _split_fixed_solute_water_blocks(solute):
    """Split an authoritative input Structure without changing atom order.

    Returns fixed non-water solute, retained input waters, and an audit. The
    retained waters are protected from generated-water ion replacement and
    later generated-water normalization.
    """
    pmd=_require_parmed()

    water_names=set(_STRUCTURAL_WATER_RESNAMES)
    original_water_residues=[
        residue for residue in solute.residues
        if residue.name.strip() in water_names
    ]
    original_water_atoms=[
        atom for residue in original_water_residues
        for atom in residue.atoms
    ]
    original_water_coordinates=(
        np.asarray(solute.coordinates,float)[
            [atom.idx for atom in original_water_atoms]
        ].copy()
        if original_water_atoms
        else np.empty((0,3),dtype=float)
    )

    water_mask=':' + ','.join(_STRUCTURAL_WATER_RESNAMES)

    # Convert AmberParm to generic Structure before splitting.
    # AmberParm.__copy__ in ParmEd 4.3.1 can fail because the
    # copied object lacks the private _ncopies attribute.
    fixed_nonwater=solute.copy(pmd.Structure)
    fixed_nonwater.strip(water_mask)

    retained_waters=solute.copy(pmd.Structure)
    retained_waters.strip('!' + water_mask)

    retained_water_residues=len(retained_waters.residues)
    retained_water_atoms=len(retained_waters.atoms)

    if retained_water_residues != len(original_water_residues):
        raise RuntimeError(
            'Structural-water split changed the retained-water residue count: '
            f'input={len(original_water_residues)}, retained={retained_water_residues}'
        )
    if retained_water_atoms != len(original_water_atoms):
        raise RuntimeError(
            'Structural-water split changed the retained-water atom count: '
            f'input={len(original_water_atoms)}, retained={retained_water_atoms}'
        )

    if original_water_atoms:
        retained_coordinates=np.asarray(
            retained_waters.coordinates,dtype=float
        )
        if not np.allclose(
            retained_coordinates,
            original_water_coordinates,
            rtol=0.0,
            atol=1.0e-8,
        ):
            raise RuntimeError(
                'Structural-water coordinates changed while splitting the input solute'
            )

    audit={
        'ordering_mode':'fixed_nonwater_then_retained_structural_water_then_environment',
        'fixed_nonwater_atoms':len(fixed_nonwater.atoms),
        'fixed_nonwater_residues':len(fixed_nonwater.residues),
        'retained_structural_water_atoms':retained_water_atoms,
        'retained_structural_water_residues':retained_water_residues,
        'retained_structural_waters_strippable':False,
        'generated_bulk_waters_strippable':True,
        'input_atom_count':len(solute.atoms),
        'split_atom_count':len(fixed_nonwater.atoms)+retained_water_atoms,
    }
    if audit['split_atom_count'] != audit['input_atom_count']:
        raise RuntimeError(
            'Fixed-solute split did not preserve the input atom count: '
            f"input={audit['input_atom_count']}, split={audit['split_atom_count']}"
        )

    return fixed_nonwater,retained_waters,audit


def _combine_protected_solute_and_environment(solute,environment,out):
    """Combine canonical blocks and prove topology/coordinate block agreement."""
    fixed_nonwater,retained_waters,audit=(
        _split_fixed_solute_water_blocks(solute)
    )
    combined=fixed_nonwater+retained_waters+environment

    n_fixed=len(fixed_nonwater.atoms)
    n_retained=len(retained_waters.atoms)
    n_environment=len(environment.atoms)

    # Structure addition must retain the exact block order. If ParmEd changes
    # it, fail here rather than write a topology/restart mismatch.
    def atom_signature(atom):
        return (
            atom.residue.name.strip(), atom.name.strip(), str(atom.type),
            round(float(atom.charge),8), round(float(atom.mass),8),
        )

    expected=(
        [atom_signature(atom) for atom in fixed_nonwater.atoms]
        + [atom_signature(atom) for atom in retained_waters.atoms]
        + [atom_signature(atom) for atom in environment.atoms]
    )
    observed=[atom_signature(atom) for atom in combined.atoms]
    if observed != expected:
        for index,(wanted,found) in enumerate(zip(expected,observed),start=1):
            if wanted != found:
                raise RuntimeError(
                    'ParmEd changed canonical component atom ordering at final '
                    f'atom {index}: expected {wanted}, observed {found}'
                )
        raise RuntimeError('ParmEd changed canonical component atom ordering')

    def coordinate_block(structure):
        if len(structure.atoms) == 0:
            return np.empty((0, 3), dtype=float)

        coordinates = np.asarray(
            structure.coordinates,
            dtype=float,
        )

        if coordinates.shape != (
            len(structure.atoms),
            3,
        ):
            raise RuntimeError(
                "Unexpected coordinate-array shape for "
                f"{len(structure.atoms)} atoms: "
                f"{coordinates.shape}"
            )

        return coordinates

    expected_coordinates=np.vstack([
        coordinate_block(fixed_nonwater),
        coordinate_block(retained_waters),
        coordinate_block(environment),
    ])
    if not np.allclose(
        np.asarray(combined.coordinates,float),
        np.asarray(expected_coordinates,float),
        rtol=0.0,
        atol=1.0e-8,
    ):
        delta=np.linalg.norm(
            np.asarray(combined.coordinates,float)
            - np.asarray(expected_coordinates,float),
            axis=1,
        )
        raise RuntimeError(
            'ParmEd changed canonical component coordinate ordering; '
            f'maximum displacement={float(delta.max()):.8f} A'
        )

    audit['blocks']={
        'fixed_nonwater':{
            'first_atom_1based':1,
            'last_atom_1based':n_fixed,
            'atom_count':n_fixed,
            'strippable':False,
        },
        'retained_structural_water':{
            'first_atom_1based':n_fixed+1 if n_retained else None,
            'last_atom_1based':n_fixed+n_retained if n_retained else None,
            'atom_count':n_retained,
            'residue_count':len(retained_waters.residues),
            'strippable':False,
        },
        'generated_environment':{
            'first_atom_1based':n_fixed+n_retained+1 if n_environment else None,
            'last_atom_1based':n_fixed+n_retained+n_environment if n_environment else None,
            'atom_count':n_environment,
            'strippable_by_subblock_metadata':True,
        },
    }
    audit['topology_coordinate_order_identical']=True

    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    (out/'structural_water_blocks.json').write_text(
        json.dumps(audit,indent=2)+'\n'
    )
    return combined,audit


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



# V0.1.4 UNIFIED POPS POPE HEADGROUP REPLACEMENT

def _convert_popc_headgroup_item(item,donors,rng,headgroup,species):
    """Replace only the POPC PC module with a PS or PE donor module."""
    molecule=item['molecule']; xyz=item['xyz']; host_atoms=molecule['atoms']
    host_pc=[(i,a) for i,a in enumerate(host_atoms) if a['resname']=='PC']
    host_map={a['name']:(i,a) for i,a in host_pc}
    donor=donors[int(rng.integers(len(donors)))]
    donor_head=[a for a in donor['atoms'] if a['resname']==headgroup]
    donor_map={a['name']:a for a in donor_head}
    alignment=['P31','O31','O32','C3','C2','C1','O21','C21','O11','C11']
    common=[name for name in alignment if name in host_map and name in donor_map]
    if len(common)<4: return None
    mobile=np.array([donor_map[name]['xyz'] for name in common])
    reference=np.array([xyz[host_map[name][0]] for name in common])
    rotation,translation=_kabsch(mobile,reference)
    retained_by_module={
        'PA': {'atoms': [], 'xyz': []},
        'OL': {'atoms': [], 'xyz': []},
    }

    for i,atom in enumerate(host_atoms):
        module=atom['resname']

        if module in retained_by_module:
            retained_by_module[module]['atoms'].append(
                dict(atom)
            )
            retained_by_module[module]['xyz'].append(
                xyz[i]
            )

    transformed_atoms=[]
    transformed_xyz=[]

    for atom in donor_head:
        copied=dict(atom)
        copied['resname']=headgroup

        transformed_atoms.append(copied)
        transformed_xyz.append(
            atom['xyz']@rotation.T+translation
        )

    # Required Lipid21 modular order:
    # PA + headgroup (PE or PS) + OL
    atoms=(
        retained_by_module['PA']['atoms']
        + transformed_atoms
        + retained_by_module['OL']['atoms']
    )

    coordinates=np.array(
        retained_by_module['PA']['xyz']
        + transformed_xyz
        + retained_by_module['OL']['xyz'],
        float,
    )
    anchors=[i for i,a in enumerate(atoms) if a['resname']==headgroup and a['name']=='P31']
    if len(anchors)!=1: return None
    result={'molecule':{'atoms':atoms,'anchor':coordinates[anchors[0]],'leaflet':item['leaflet']},
            'xyz':coordinates,'anchor':coordinates[anchors[0]],
            'elements':[a['element'] for a in atoms], 'leaflet':item['leaflet'],'species':species}
    if 'canonical_site_id' in item: result['canonical_site_id']=item['canonical_site_id']
    return result


def _place_exact_headgroup_target(items,eligible,target,donors,fixed_cells,cfg,rng,leaf,headgroup,species):
    """Unified exact POPC-to-POPS/POPE replacement with final periodic checks."""
    candidates=[index for index in eligible if items[index]['species']=='POPC']
    rng.shuffle(candidates); placed=[]; failures=[]
    donor_pool=[m for m in donors if m['leaflet']==leaf] or donors
    key=species.lower()+'_donors_per_site'
    donor_limit=min(len(donor_pool),int(cfg.get('complex_membrane',{}).get(key,64)))
    for index in candidates:
        if len(placed)>=target: break
        order=list(rng.permutation(len(donor_pool)))[:donor_limit]; best=None
        for donor_index in order:
            attempt=_convert_popc_headgroup_item(items[index],[donor_pool[donor_index]],rng,headgroup,species)
            if attempt is None: continue
            valid,score,minima=_validate_replacement_item(attempt,index,items,fixed_cells,cfg)
            if not valid: continue
            if best is None or score<best[0]: best=(score,attempt,minima)
            if score<=1e-12: break
        if best is None:
            failures.append({'site_index':index,'reason':f'no_valid_{species.lower()}_conversion'}); continue
        score,attempt,minima=best; items[index]=attempt
        placed.append({'site_index':index,'score':float(score),'minimum_distances':minima})
    if len(placed)!=target:
        raise RuntimeError(f'Exact composition failed: requested {target} {species} in {leaf}, placed {len(placed)} after testing {len(candidates)} POPC sites.')
    return placed,failures


def _allocate_leaflet_counts(total,fractions,order=('POPC','POPE','POPS','CHL')):
    """Largest-remainder integer allocation with a deterministic tie break."""
    clean={name:float(fractions.get(name,0.0)) for name in order}
    if any(v<0 for v in clean.values()): raise ValueError(f'Negative leaflet fraction: {clean}')
    s=sum(clean.values())
    if abs(s-1.0)>1e-8: raise ValueError(f'Leaflet fractions must sum to 1.0, got {s}: {clean}')
    raw={name:total*clean[name] for name in order}; counts={name:int(np.floor(raw[name])) for name in order}
    remaining=total-sum(counts.values())
    ranked=sorted(
        order,
        key=lambda name:(
            -round(raw[name] - counts[name], 12),
            order.index(name),
        ),
    )
    for name in ranked[:remaining]: counts[name]+=1
    return counts


def _composition_fraction(settings,cfg,leaf,species,default=0.0):
    modern=cfg.get('membrane',{}).get('composition',{}).get(leaf,{})
    if species in modern: return float(modern[species])
    old=settings.get('leaflets',{}).get(leaf,{})
    return float(old.get(species+'_fraction',settings.get(species+'_fraction',default)))


def _packaged_lipid_reference(species):
    return Path(__file__).resolve().parent/'data'/'reference_library'/species/'reference.pdb'

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
    """Validate a replacement against every other molecule with periodic MIC."""
    environment=[item for index,item in enumerate(items) if index!=host_index]
    box=cfg.get('_periodic_box_angstrom')
    if box is not None:
        return _periodic_candidate_valid(candidate['xyz'],candidate['elements'],environment,box,cfg)
    lipid_cells=_build_atom_cells(environment,4.0)
    return _candidate_valid(
        candidate['xyz'],candidate['elements'],fixed_cells,lipid_cells,4.0,
        float(cfg.get('packing',{}).get('fixed_heavy_lipid_heavy_angstrom',2.0)),
        float(cfg.get('packing',{}).get('fixed_heavy_lipid_hydrogen_angstrom',1.2)),
        1.5,1.05,float(cfg.get('packing',{}).get('lipid_hydrogen_hydrogen_angstrom',1.15)))


# V0.1.5 GLOBAL EXACT CHOLESTEROL PLACEMENT

def _cholesterol_periodic_box_xy(cfg):
    """Resolve the commensurate membrane XY dimensions for site ordering."""
    requested=(cfg.get('membrane_only',{}).get('box_xy_angstrom')
               or cfg.get('membrane',{}).get('box_xy_angstrom'))
    if not requested: return None
    reference=_packaged_lipid_reference('POPC')
    line=next((line for line in Path(reference).read_text(errors='replace').splitlines()
               if line.startswith('CRYST1')),None)
    if line is None: return np.asarray(requested[:2],float)
    reference_xy=np.array([float(line[6:15]),float(line[15:24])],float)
    repeats=np.maximum(1,np.rint(np.asarray(requested[:2],float)/reference_xy).astype(int))
    return reference_xy*repeats


def _cholesterol_spread_candidate_order(candidates,items,cfg,rng):
    """Seeded periodic farthest-point order for dense CHL replacement.

    The order is spatially distributed, so early accepted cholesterol sites do
    not cluster and block the remaining exact target. All eligible sites remain
    in the returned order, preserving exact-composition semantics.
    """
    candidates=list(candidates)
    if len(candidates)<2: return candidates
    points=np.array([items[index]['anchor'][:2] for index in candidates],float)
    box=_cholesterol_periodic_box_xy(cfg)
    chosen=[int(rng.integers(len(candidates)))]
    remaining=set(range(len(candidates))); remaining.remove(chosen[0])
    min_distance=np.full(len(candidates),np.inf,float)
    while remaining:
        last=chosen[-1]
        delta=points-points[last]
        if box is not None:
            delta-=box*np.rint(delta/box)
        distance2=np.einsum('ij,ij->i',delta,delta)
        min_distance=np.minimum(min_distance,distance2)
        maximum=max(min_distance[index] for index in remaining)
        tied=[index for index in remaining if abs(min_distance[index]-maximum)<=1e-10]
        next_index=tied[int(rng.integers(len(tied)))]
        chosen.append(next_index); remaining.remove(next_index)
    return [candidates[index] for index in chosen]


def _place_exact_cholesterol_target_single_attempt(items,eligible,target,donors,fixed_cells,cfg,rng,leaf):
    """Try seed-ordered POPC hosts until the exact CHL target is committed."""
    candidates=list(eligible); candidates=_cholesterol_spread_candidate_order(candidates,items,cfg,rng)
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
        candidate['canonical_site_id']=items[index].get('canonical_site_id',(index,))
        items[index]=candidate
        placed.append({'site_index':index,'score':float(score),'minimum_distances':minima})
    if len(placed)!=target:
        raise RuntimeError(
            f'Exact composition failed: requested {target} CHL in {leaf} leaflet, '
            f'placed {len(placed)} after testing {len(candidates)} eligible POPC sites.'
        )
    return placed,failures



def _place_exact_cholesterol_target(items,eligible,target,donors,fixed_cells,cfg,rng,leaf):
    """Place an exact CHL target with distributed sites and pristine retries."""
    if target<=0: return [],[]
    settings=cfg.get('complex_membrane',{})
    attempts=max(1,int(settings.get('cholesterol_global_attempts',12)))
    original=copy.deepcopy(items)
    errors=[]
    for attempt_number in range(1,attempts+1):
        items[:]=copy.deepcopy(original)
        attempt_seed=int(rng.integers(0,np.iinfo(np.uint32).max,dtype=np.uint32))
        attempt_rng=np.random.default_rng(attempt_seed)
        try:
            placed,failures=_place_exact_cholesterol_target_single_attempt(
                items,eligible,target,donors,fixed_cells,cfg,attempt_rng,leaf)
            return placed,failures
        except RuntimeError as error:
            errors.append({'attempt':attempt_number,'seed':attempt_seed,'error':str(error)})
    items[:]=original
    summary=' | '.join(f"attempt {item['attempt']}: {item['error']}" for item in errors[-3:])
    raise RuntimeError(
        f'Exact cholesterol composition failed after {attempts} distributed attempts '
        f'for {leaf} leaflet (target {target}). Recent failures: {summary}'
    )

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
        attempt['canonical_site_id']=items[index].get('canonical_site_id',(index,))
        items[index]=attempt
        placed.append({'site_index':index,'score':float(score),'minimum_distances':minima})
    if len(placed)!=target:
        raise RuntimeError(
            f'Exact composition failed: requested {target} POPS in {leaf} leaflet, '
            f'placed {len(placed)} after testing {len(candidates)} eligible POPC sites.'
        )
    return placed,failures


def _apply_complex_composition(items,fixed,lower,upper,cfg,rng):
    """Apply exact CHL/custom/POPS/POPE composition with shared headgroup logic."""
    settings=cfg.get('complex_membrane',{})
    modern=cfg.get('membrane',{}).get('composition',{})
    if not settings.get('enabled',False) and not modern: return items,{}
    membrane_fixed=[atom for atom in fixed if lower-8<=atom['xyz'][2]<=upper+8]
    shell=float(settings.get('protein_exclusion_shell_angstrom',8.0))
    chl_needed=any(_composition_fraction(settings,cfg,leaf,'CHL',0.0)>0 for leaf in ('upper','lower'))
    chl_ref=_parse_chl_reference(settings.get('CHL_reference',_packaged_lipid_reference('CHL').with_name('conformers.pdb'))) if chl_needed else []
    custom_settings=_custom_sterol_settings(cfg); custom_template=_custom_sterol_template(custom_settings) if custom_settings else None
    pops_needed=any(_composition_fraction(settings,cfg,leaf,'POPS',0.0)>0 for leaf in ('upper','lower'))
    pope_needed=any(_composition_fraction(settings,cfg,leaf,'POPE',0.0)>0 for leaf in ('upper','lower'))
    pops_ref=_parse_modular_lipid_reference(settings.get('POPS_reference',_packaged_lipid_reference('POPS')),'PS','POPS') if pops_needed else []
    pope_ref=_parse_modular_lipid_reference(settings.get('POPE_reference',_packaged_lipid_reference('POPE')),'PE','POPE') if pope_needed else []
    fixed_cells=build_fixed_cells(fixed,4.0)
    report={'order':['CHL','CUSTOM_STEROL','POPS','POPE'],'protected_shell_angstrom':shell,'leaflets':{}}
    for leaf in ('upper','lower'):
        idx=[i for i,item in enumerate(items) if item['leaflet']==leaf]; total=len(idx)
        fractions={name:_composition_fraction(settings,cfg,leaf,name,1.0 if name=='POPC' else 0.0) for name in ('POPC','POPE','POPS','CHL')}
        custom_fraction=0.0
        if custom_settings:
            custom_fraction=float(custom_settings.get('leaflets',{}).get(leaf,{}).get('fraction',custom_settings.get('fraction',0.0)))
        # Legacy configurations omit POPC; infer it after fixed species are read.
        modern_leaf=modern.get(leaf,{})
        old_leaf=settings.get('leaflets',{}).get(leaf,{})
        if 'POPC' not in modern_leaf and 'POPC_fraction' not in old_leaf:
            fractions['POPC']=1.0-fractions['POPE']-fractions['POPS']-fractions['CHL']-custom_fraction
        alloc_fractions=dict(fractions)
        if custom_fraction:
            alloc_fractions['CUSTOM_STEROL']=custom_fraction
            order=('POPC','POPE','POPS','CHL','CUSTOM_STEROL')
        else: order=('POPC','POPE','POPS','CHL')
        counts=_allocate_leaflet_counts(total,alloc_fractions,order)
        eligible=[i for i in idx if items[i]['species']=='POPC' and _distance_to_fixed_xy(items[i]['anchor'],membrane_fixed)>=shell]
        replacements=counts.get('CHL',0)+counts.get('CUSTOM_STEROL',0)+counts.get('POPS',0)+counts.get('POPE',0)
        if len(eligible)<replacements: raise RuntimeError(f'Exact composition impossible in {leaf}: {len(eligible)} eligible sites for {replacements} replacements')
        chl,chl_bad=(_place_exact_cholesterol_target(items,eligible,counts['CHL'],chl_ref,fixed_cells,cfg,rng,leaf) if counts['CHL'] else ([],[]))
        custom,custom_bad=(_place_exact_custom_sterol_target(items,eligible,counts['CUSTOM_STEROL'],custom_template,fixed_cells,cfg,rng,leaf,custom_settings) if counts.get('CUSTOM_STEROL',0) else ([],[]))
        pops,pops_bad=(_place_exact_headgroup_target(items,eligible,counts['POPS'],pops_ref,fixed_cells,cfg,rng,leaf,'PS','POPS') if counts['POPS'] else ([],[]))
        pope,pope_bad=(_place_exact_headgroup_target(items,eligible,counts['POPE'],pope_ref,fixed_cells,cfg,rng,leaf,'PE','POPE') if counts['POPE'] else ([],[]))
        species=['POPC','POPE','POPS','CHL']+([custom_template['species']] if custom_template else [])
        actual={name:sum(items[i]['species']==name for i in idx) for name in species}
        expected={name:counts.get(name,0) for name in ('POPC','POPE','POPS','CHL')}
        if custom_template: expected[custom_template['species']]=counts.get('CUSTOM_STEROL',0)
        if actual!=expected: raise RuntimeError(f'Composition verification failed in {leaf}: expected {expected}, got {actual}')
        report['leaflets'][leaf]={'total':total,'requested_fractions':alloc_fractions,'target':expected,'actual':actual,'eligible_bulk_sites':len(eligible),
          'cholesterol_placements':chl,'cholesterol_failed_sites':chl_bad,'custom_sterol_placements':custom,'custom_sterol_failed_sites':custom_bad,
          'pops_conversions':pops,'pops_failed_sites':pops_bad,'pope_conversions':pope,'pope_failed_sites':pope_bad}
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


def _water_overlaps(xyz,elements,cells,oxygen_cut,hydrogen_cut,virtual_cut,cell_width=4.0,box=None):
    if box is not None:
        others=[other for bucket in cells.values() for other,_ in bucket]
        if not others: return False
        others=np.asarray(others,float); box=np.asarray(box,float)
        for point,element in zip(xyz,elements):
            cutoff=oxygen_cut if element=='O' else (hydrogen_cut if element=='H' else virtual_cut)
            delta=_mic_delta(others-point,box)
            if np.any(np.einsum('ij,ij->i',delta,delta)<cutoff*cutoff): return True
        return False
    for point,element in zip(xyz,elements):
        cutoff=oxygen_cut if element=='O' else (hydrogen_cut if element=='H' else virtual_cut)
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
                        oxygen_cut,hydrogen_cut,virtual_cut,box=box
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
                    # V0.16.0 register accepted OPC water in collision cells
                    for accepted_point,accepted_element in zip(xyz,elements):
                        accepted_key=tuple(np.floor(accepted_point/4.0).astype(int))
                        cells[accepted_key].append((accepted_point,accepted_element))

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
            periodic_box=cfg.get('_periodic_box_angstrom')
            if fixed_points.size:
                delta=fixed_points-point
                if periodic_box is not None: delta=_mic_delta(delta,np.asarray(periodic_box,float))
                if np.min(np.linalg.norm(delta,axis=1))<minimum_fixed: continue
            if selected:
                delta=np.asarray([other[1] for other in selected])-point
                if periodic_box is not None: delta=_mic_delta(delta,np.asarray(periodic_box,float))
                if np.min(np.linalg.norm(delta,axis=1))<minimum_ion: continue
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
    solute,structural_water_cleanup,prepared_solute_top,prepared_solute_rst=(
        _prepare_authoritative_solute_pair(solute_prmtop,solute_rst7,cfg,out)
    )
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
    # Reload the authoritative input pair immediately before final assembly.
    # The long-lived working solute can acquire stale ParmEd molecule/cache
    # state during packing and topology preparation. Never use it for writing
    # the final topology/restart pair.
    fresh_solute=pmd.load_file(
        str(prepared_solute_top),xyz=str(prepared_solute_rst)
    )
    final_structural_water_cleanup=structural_water_cleanup
    combined,structural_water_audit=(
        _combine_protected_solute_and_environment(
            fresh_solute,final_environment,out
        )
    )
    combined.box=np.array([box[0],box[1],box[2],90.0,90.0,90.0])
    solute_preservation_audit=_verify_solute_prefix(fresh_solute,combined)
    final_charge=_rounded_integer_charge(sum(float(a.charge) for a in combined.atoms),tolerance)
    if ion_cfg.get('enabled',False) and final_charge!=0:
        raise RuntimeError(f'Final topology charge is {final_charge} e, not zero')
    topology_name=top_cfg.get('output_prmtop','topology.prmtop')
    coordinate_name=top_cfg.get('output_rst7','coordinates.rst7')
    expected_final_atoms=len(combined.atoms)

    # ParmEd's direct periodic Structure writer can attempt an implicit
    # molecule reorder and corrupt the topology atom count. First write a
    # stable nonperiodic pair from independent generic Structure copies.
    temporary_topology=out/'final_unboxed_intermediate.prmtop'
    temporary_coordinates=out/'final_unboxed_intermediate.rst7'

    for temporary_path in (temporary_topology,temporary_coordinates):
        temporary_path.unlink(missing_ok=True)

    unboxed_topology=combined.copy(pmd.Structure)
    unboxed_coordinates=combined.copy(pmd.Structure)
    unboxed_topology.box=None
    unboxed_coordinates.box=None

    with warnings.catch_warnings(record=True) as unboxed_warnings:
        warnings.simplefilter('always')
        unboxed_topology.save(
            str(temporary_topology),
            format='amber',
            overwrite=True,
        )

    molecule_warnings=[
        str(item.message)
        for item in unboxed_warnings
        if (
            'not contiguous' in str(item.message).lower()
            or 'reorder' in str(item.message).lower()
        )
    ]
    if molecule_warnings:
        raise RuntimeError(
            'Ordered final structure triggered molecular reordering during '
            'the nonperiodic topology write: '
            + ' | '.join(molecule_warnings)
        )

    unboxed_coordinates.save(
        str(temporary_coordinates),
        format='rst7',
        overwrite=True,
    )

    temporary_written_topology=pmd.load_file(str(temporary_topology))
    temporary_restart_atoms=int(
        temporary_coordinates.read_text().splitlines()[1].split()[0]
    )
    if (
        len(temporary_written_topology.atoms)!=expected_final_atoms
        or temporary_restart_atoms!=expected_final_atoms
    ):
        raise RuntimeError(
            'Intermediate topology/restart atom-count mutation: '
            f'in-memory={expected_final_atoms}, '
            f'topology={len(temporary_written_topology.atoms)}, '
            f'restart={temporary_restart_atoms}'
        )

    # Reload as AmberParm, explicitly construct Amber molecule metadata,
    # assign the periodic box, and write the final exact pair from this same
    # ordered object. Do not call AmberParm.copy(); ParmEd 4.3.1 can create
    # an incomplete copy without parm_data["POINTERS"].
    periodic=pmd.load_file(
        str(temporary_topology),
        xyz=str(temporary_coordinates),
    )

    from parmed.tools.actions import setMolecules

    with warnings.catch_warnings(record=True) as molecule_metadata_warnings:
        warnings.simplefilter('always')
        setMolecules(periodic).execute()

    metadata_reorder_warnings=[
        str(item.message)
        for item in molecule_metadata_warnings
        if (
            'not contiguous' in str(item.message).lower()
            or 'reorder' in str(item.message).lower()
        )
    ]
    if metadata_reorder_warnings:
        raise RuntimeError(
            'setMolecules attempted to reorder the normalized final system: '
            + ' | '.join(metadata_reorder_warnings)
        )

    periodic.box=np.array(
        [box[0],box[1],box[2],90.0,90.0,90.0],
        dtype=float,
    )

    with warnings.catch_warnings(record=True) as periodic_warnings:
        warnings.simplefilter('always')
        periodic.save(
            str(out/topology_name),
            overwrite=True,
        )

    periodic_reorder_warnings=[
        str(item.message)
        for item in periodic_warnings
        if (
            'not contiguous' in str(item.message).lower()
            or 'reorder' in str(item.message).lower()
        )
    ]
    if periodic_reorder_warnings:
        raise RuntimeError(
            'Final periodic topology write attempted molecular reordering: '
            + ' | '.join(periodic_reorder_warnings)
        )

    periodic.save(
        str(out/coordinate_name),
        format='rst7',
        overwrite=True,
    )

    written_topology=pmd.load_file(str(out/topology_name))
    restart_second_line=(out/coordinate_name).read_text().splitlines()[1]
    written_restart_atoms=int(restart_second_line.split()[0])
    if len(written_topology.atoms)!=expected_final_atoms or written_restart_atoms!=expected_final_atoms:
        raise RuntimeError(
            'Final topology/restart atom-count mutation during write: '
            f'in-memory={expected_final_atoms}, topology={len(written_topology.atoms)}, '
            f'restart={written_restart_atoms}'
        )

    exact_final_system=pmd.load_file(
        str(out/topology_name),
        xyz=str(out/coordinate_name),
    )
    if len(exact_final_system.atoms)!=expected_final_atoms:
        raise RuntimeError(
            'Exact final topology/restart reload changed the atom count: '
            f'expected={expected_final_atoms}, '
            f'reloaded={len(exact_final_system.atoms)}'
        )

    expected_box=np.array(
        [box[0],box[1],box[2],90.0,90.0,90.0],
        dtype=float,
    )
    if exact_final_system.box is None:
        raise RuntimeError(
            'Exact final topology/restart pair has no periodic box'
        )
    if not np.allclose(
        np.asarray(exact_final_system.box,dtype=float),
        expected_box,
        rtol=0.0,
        atol=1.0e-3,
    ):
        raise RuntimeError(
            'Exact final topology/restart box differs from the builder box: '
            f'expected={expected_box.tolist()}, '
            f'found={np.asarray(exact_final_system.box,dtype=float).tolist()}'
        )

    validation_cfg=cfg.get('final_validation',{})
    if (
        structural_water_cleanup.get('hard_clashing_structural_waters',0)>0
        and structural_water_cleanup.get('policy')=='keep'
    ):
        raise RuntimeError(
            'Structural waters with hard clashes were kept for diagnosis; '
            'SIMULATION_READY is not permitted'
        )

    final_validation_report=validate_final_system(
        out/topology_name,out/coordinate_name,
        unrelated_cutoff_angstrom=float(validation_cfg.get('unrelated_atom_cutoff_angstrom',0.5)),
        heavy_cutoff_angstrom=float(validation_cfg.get('nonbonded_heavy_cutoff_angstrom',1.0)),
        bond_offset_angstrom=float(validation_cfg.get('bond_length_offset_angstrom',0.5)),
        energy_threshold_kj_mol=float(validation_cfg.get('energy_threshold_kj_mol',1.0e12)),
        cpptraj_cutoff_angstrom=float(validation_cfg.get('cpptraj_cutoff_angstrom',10.0)),
        report_limit=int(validation_cfg.get('report_limit',100)),
        report_path=out/'final_validation.json')
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
            'solute_bonds_preserved':solute_preservation_audit['bond_definitions_preserved'],
            'final_validation':final_validation_report}
    return final_items,result,combined

# V0.15.0 MEMBRANE-ONLY BUILDER
def _build_membrane_only_host(molecules,reference_box,requested_box_xy,rng):
    """Build the nearest complete periodic reference supercell.

    The requested X/Y dimensions are preferences. Each final dimension is an
    integer multiple of the validated reference CRYST1 vector. Every molecule
    in the selected supercell is retained, so no vacancies or incompatible
    periodic crop are introduced.
    """
    requested=np.asarray(requested_box_xy,float)
    reference_box=np.asarray(reference_box,float)
    if requested.shape!=(2,) or np.any(requested<=0):
        raise ValueError('membrane_only.box_xy_angstrom must contain two positive values')
    if reference_box.size<2 or np.any(reference_box[:2]<=0):
        raise ValueError('Validated POPC reference has invalid CRYST1 dimensions')
    repeats=np.maximum(1,np.rint(requested/reference_box[:2]).astype(int))
    # If nearest rounding undershoots by more than half a cell, ceil naturally
    # wins. rint provides the mathematically nearest commensurate supercell.
    actual=reference_box[:2]*repeats
    selected=[]
    for ix in range(int(repeats[0])):
        for iy in range(int(repeats[1])):
            shift=np.array([ix*reference_box[0],iy*reference_box[1],0.0])
            for molecule_index,molecule in enumerate(molecules):
                xyz,anchor,elements=_atom_arrays_from_molecule(molecule,shift=shift)
                selected.append({'molecule':molecule,'xyz':xyz,'anchor':anchor,
                    'elements':elements,'leaflet':molecule['leaflet'],'species':'POPC',
                    'canonical_site_id':(int(ix),int(iy),int(molecule_index))})
    counts=Counter(item['leaflet'] for item in selected)
    if counts['upper']!=counts['lower']:
        raise RuntimeError(f'Periodic host is leaflet-asymmetric: {dict(counts)}')
    # Shift the complete supercell so X/Y are centered near zero. This is a
    # uniform translation and preserves the validated periodic seam.
    shift_xy=-0.5*actual
    for item in selected:
        item['xyz'][:,:2]+=shift_xy
        item['anchor'][:2]+=shift_xy
    report={
        'box_mode':'nearest_commensurate_reference_supercell',
        'requested_box_xy_angstrom':[float(v) for v in requested],
        'reference_box_xy_angstrom':[float(v) for v in reference_box[:2]],
        'reference_repeats_xy':[int(v) for v in repeats],
        'actual_box_xy_angstrom':[float(v) for v in actual],
        'adjustment_xy_angstrom':[float(v) for v in actual-requested],
        'selected_sites_per_leaflet':int(counts['upper']),
        'selected_total_sites':len(selected),
    }
    return selected,actual,report


def _mic_delta(delta,box):
    box=np.asarray(box,float)
    return delta-box*np.rint(delta/box)


def _periodic_candidate_valid(candidate_xyz,candidate_elements,environment,box,cfg):
    """Vectorized exact minimum-image validation against all other molecules."""
    box=np.asarray(box,float)
    other_xyz=[]; other_elements=[]
    for item in environment:
        other_xyz.extend(np.asarray(item['xyz'],float))
        other_elements.extend(item.get('elements',[atom['element'] for atom in item['molecule']['atoms']]))
    if not other_xyz:
        return True,0.0,{'lipid_heavy_heavy':999.0,'lipid_heavy_hydrogen':999.0,'lipid_hydrogen_hydrogen':999.0}
    other_xyz=np.asarray(other_xyz,float); other_elements=np.asarray(other_elements)
    hh=1.5; hx=1.05; xx=float(cfg.get('packing',{}).get('lipid_hydrogen_hydrogen_angstrom',1.15))
    minima={'fixed_heavy_lipid_heavy':999.0,'fixed_heavy_lipid_hydrogen':999.0,
            'lipid_heavy_heavy':999.0,'lipid_heavy_hydrogen':999.0,'lipid_hydrogen_hydrogen':999.0}
    score=0.0
    for point,element in zip(np.asarray(candidate_xyz,float),candidate_elements):
        delta=_mic_delta(other_xyz-point,box)
        distances=np.linalg.norm(delta,axis=1)
        if element=='H':
            hmask=other_elements=='H'; xmask=~hmask
            if np.any(hmask): minima['lipid_hydrogen_hydrogen']=min(minima['lipid_hydrogen_hydrogen'],float(distances[hmask].min()))
            if np.any(xmask): minima['lipid_heavy_hydrogen']=min(minima['lipid_heavy_hydrogen'],float(distances[xmask].min()))
            if np.any(hmask & (distances<xx)) or np.any(xmask & (distances<hx)): return False,1e9,minima
        else:
            hmask=other_elements=='H'; xmask=~hmask
            if np.any(hmask): minima['lipid_heavy_hydrogen']=min(minima['lipid_heavy_hydrogen'],float(distances[hmask].min()))
            if np.any(xmask): minima['lipid_heavy_heavy']=min(minima['lipid_heavy_heavy'],float(distances[xmask].min()))
            if np.any(hmask & (distances<hx)) or np.any(xmask & (distances<hh)): return False,1e9,minima
    return True,score,minima

def _canonicalize_environment_item_order(items):
    """Return deterministic molecular order independent of seeded site choices.

    Membrane molecules are grouped by topology identity, then leaflet, then
    stable host-site ID. Solvent remains after the membrane in deterministic
    generation order. Ions are appended after replacing waters.
    """
    species_rank={'POPC':0,'POPS':1,'POPE':2,'CHL':3}
    leaflet_rank={'upper':0,'lower':1,'positive':0,'negative':1}
    membrane=[]; solvent=[]; other=[]
    for generation_index,item in enumerate(items):
        species=str(item.get('species',''))
        if species=='OPC':
            solvent.append((generation_index,item))
        elif species in {'K+','Na+','Cl-','K','NA','CL'}:
            other.append((generation_index,item))
        else:
            membrane.append((generation_index,item))

    def membrane_key(entry):
        generation_index,item=entry
        species=str(item.get('species',''))
        rank=species_rank.get(species,10)
        # Custom residues are ordered deterministically by residue/species name.
        custom_name=species if rank==10 else ''
        leaflet=leaflet_rank.get(str(item.get('leaflet','')),9)
        site=tuple(item.get('canonical_site_id',(generation_index,)))
        return (rank,custom_name,leaflet,site)

    membrane.sort(key=membrane_key)
    solvent.sort(key=lambda entry:entry[0])
    other.sort(key=lambda entry:(str(entry[1].get('species','')),entry[0]))
    ordered=[item for _,item in membrane+solvent+other]
    return ordered,{
        'mode':'species_leaflet_canonical_site',
        'membrane_items':len(membrane),
        'solvent_items':len(solvent),
        'other_items':len(other),
        'species_order':['POPC','POPS','POPE','CHL','custom_by_species','OPC','ions'],
        'seed_independent_topology_order':True,
    }


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
    validation_cfg=cfg.get('final_validation',{})
    final_validation_report=validate_final_system(
        out/topology_name,out/coordinate_name,
        unrelated_cutoff_angstrom=float(validation_cfg.get('unrelated_atom_cutoff_angstrom',0.5)),
        heavy_cutoff_angstrom=float(validation_cfg.get('nonbonded_heavy_cutoff_angstrom',1.0)),
        bond_offset_angstrom=float(validation_cfg.get('bond_length_offset_angstrom',0.5)),
        energy_threshold_kj_mol=float(validation_cfg.get('energy_threshold_kj_mol',1.0e12)),
        cpptraj_cutoff_angstrom=float(validation_cfg.get('cpptraj_cutoff_angstrom',10.0)),
        report_limit=int(validation_cfg.get('report_limit',100)),
        report_path=out/'final_validation.json')
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
        'final_validation':final_validation_report,
    },final_system


def _build_membrane_only(config_path,cfg):
    settings=cfg.get('membrane_only',{})
    box_xy=np.asarray(settings.get('box_xy_angstrom',[110.0,110.0]),float)
    # Requested counts are advisory in commensurate mode. The complete
    # reference supercell determines the simulation-safe exact count.
    rng=np.random.default_rng(int(cfg.get('build',{}).get('seed',1001)))
    reference=Path(cfg['reference_library']['POPC']).expanduser().resolve()
    molecules,reference_box=parse_reference(reference)
    assembled,actual_box_xy,host_report=_build_membrane_only_host(
        molecules,reference_box,box_xy,rng)
    target_per_leaflet=int(host_report['selected_sites_per_leaflet'])
    target_total=int(host_report['selected_total_sites'])
    lower=-float(settings.get('half_thickness_angstrom',17.9))
    upper=+float(settings.get('half_thickness_angstrom',17.9))
    box=np.array([actual_box_xy[0],actual_box_xy[1],float(cfg.get('membrane',{}).get('box_z_angstrom',96.0))])
    fixed=[]
    cfg['_periodic_box_angstrom']=[float(box[0]),float(box[1]),float(box[2])]
    assembled,composition_report=_apply_complex_composition(
        assembled,fixed,lower,upper,cfg,rng)
    assembled,canonical_order_report=_canonicalize_environment_item_order(assembled)
    assembled,box,solvent_report=_apply_opc_solvation(
        assembled,fixed,box,cfg)
    out=Path(cfg['output']['directory']); out.mkdir(parents=True,exist_ok=True)
    assembled,topology_report,final_system=_prepare_membrane_only_topology_and_ions(
        assembled,box,cfg,out,rng)
    output_pdb=out/cfg.get('output',{}).get('pdb_name','system_membrane_solvated.pdb')
    final_system.save(str(output_pdb),overwrite=True)
    counts=Counter((item['species'],item['leaflet']) for item in assembled)
    report={
        'version':'0.16.0-cpptraj-final-gate','build_mode':'membrane_only',
        'seed':int(cfg.get('build',{}).get('seed',1001)),
        'box_angstrom':[float(value) for value in box],
        'host_selection':host_report,
        'canonical_ordering':canonical_order_report,
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
    assembled,canonical_order_report=_canonicalize_environment_item_order(assembled)
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
      'version':'0.16.0-cpptraj-final-gate','input_pdb':str(input_path),'input_sha256':sha256(input_path),
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
      'canonical_ordering':canonical_order_report,
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
