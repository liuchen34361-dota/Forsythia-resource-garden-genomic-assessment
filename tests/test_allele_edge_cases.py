#!/usr/bin/env python3
"""Independent hand-calculated checks for SNP states and partial GT."""
import csv,gzip,json,shutil,subprocess,sys,tempfile,importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'scripts'))
from analysis_utils import connected_components
assert connected_components(['x10','x2','x3','single'], [('x10','x2'),('x2','x3')]) == [['single'],['x10','x2','x3']]
spec=importlib.util.spec_from_file_location('portfolio',P/'scripts/10_final_wild_portfolios.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert m.called_gt('1/.')==[1] and m.called_gt('.|0')==[0]
try:m.called_gt('0/2')
except ValueError:pass
else:raise AssertionError('Out-of-range GT accepted')
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp)/'project';shutil.copytree(P/'tests/synthetic_project',root)
 cfg=root/'config/analysis_config.json';c=json.loads(cfg.read_text());c['project_root']=str(root);cfg.write_text(json.dumps(c))
 vcf=root/'01_qc/04.samples_filtered.PASS.snp.vcf.gz';s=gzip.open(vcf,'rt').read()
 # Current E1/E2/E3/E4; E2 removed. Star only in removed E2 must not count.
 def row(pos,ref,alt,rg):return '\t'.join(['chr1',str(pos),'.',ref,alt,'.','PASS','.','GT']+['0/0']*6+rg)+'\n'
 s+=row(7,'T','C,*',['1/.','2/2','0/0','0/0'])
 s+=row(8,'a','g',['0/0','1/.','0/0','0/0'])
 s+=row(9,'A','A',['1/1','1/1','0/0','0/0'])
 with gzip.open(vcf,'wt') as f:f.write(s)
 for script in ['05_PIHAT_components_and_core.py','06_core_SNP_ALT_coverage.py']:
  subprocess.run([sys.executable,str(P/'scripts'/script),'--config',str(cfg)],check=True,stdout=subprocess.DEVNULL)
 with (root/'results/06_core_SNP_ALT_coverage/01.deredundancy_ALT_coverage_summary.tsv').open() as f:r=next(csv.DictReader(f,delimiter='\t'))
 for k,v in {'Current_ex_situ_observed_ALT_states':6,'ALT_states_retained_in_core':4,'ALT_states_lost_after_deredundancy':2,'Current_ex_situ_ALT_copies':11,'Current_ALT_copies_covered_by_core_presence':9}.items():assert int(r[k])==v,(k,r[k],v)
 again=subprocess.run([sys.executable,str(P/'scripts/06_core_SNP_ALT_coverage.py'),'--config',str(cfg)],capture_output=True,text=True)
 assert again.returncode!=0 and 'Output step is not empty' in again.stderr
print('[PASS] Mixed SNP/star, partial GT, lowercase bases, identical REF/ALT exclusion, transitive components and overwrite protection.')
