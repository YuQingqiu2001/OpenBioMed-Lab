"""从本地原始研究目录复现权重整合。

公开包不附带原始教师与 H&E 特征，不自动下载。
输出工作目录必须位于只读原始研究目录之外。
"""
from pathlib import Path
import argparse,os,runpy,sys

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--study-root',type=Path,required=True,help='冻结研究目录，含 references/、checkpoints/checkpoint_registry.csv 和 runs/')
 p.add_argument('--work-root',type=Path,required=True,help='独立整合工作目录，保存目标缓存与学生候选')
 p.add_argument('--phase',choices=['cache','targets','student'],required=True,help='cache 整理输入，targets 集成目标，student 训练与选择')
 a,remaining=p.parse_known_args()
 study=a.study_root.resolve();work=a.work_root.resolve()
 if study==work or study in work.parents:raise ValueError('整合输出必须位于原始研究目录之外')
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
