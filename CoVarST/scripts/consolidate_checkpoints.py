"""Optional reproducible consolidation from a local original study tree.

Original teachers and original H&E features are not distributed. Nothing is
downloaded. The work directory must be separate from the read-only study root.
"""
from pathlib import Path
import argparse,os,runpy,sys

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--study-root',type=Path,required=True,help='Frozen study tree with references/, checkpoints/checkpoint_registry.csv and runs/')
 p.add_argument('--work-root',type=Path,required=True,help='Fresh local consolidation workspace; holds target caches and student candidates')
 p.add_argument('--phase',choices=['cache','targets','student'],required=True)
 a,remaining=p.parse_known_args()
 study=a.study_root.resolve();work=a.work_root.resolve()
 if study==work or study in work.parents:raise ValueError('Keep consolidation outputs outside the original scientific study tree')
 os.environ['COVARST_STUDY_ROOT']=str(study);os.environ['COVARST_BUILD_ROOT']=str(work)
 root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
 from covarst.runtime import activate,ENGINE
 activate('he','consolidation')
 script={'cache':'prepare_distillation.py','targets':'ensemble_targets.py','student':'train_ensemble_student.py'}[a.phase]
 path=ENGINE/'consolidation'/script
 if remaining[:1]==['--']:remaining.pop(0)
 sys.argv=[str(path)]+remaining
 runpy.run_path(str(path),run_name='__main__')
if __name__=='__main__':main()
