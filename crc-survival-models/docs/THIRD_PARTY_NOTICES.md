# External dependencies / 外部依赖

This repository does not contain or redistribute UNI or CellViT++ source code,
checkpoints, classifier heads, or derived large backbone files.

本仓库不包含也不再分发 UNI 或 CellViT++ 的源码、主干权重、分类头或衍生大模型文件。

## Official upstream repositories / 官方原库

- CellViT++: https://github.com/TIO-IKIM/CellViT-plus-plus
- UNI: https://github.com/mahmoodlab/UNI

Users must obtain any required third-party files directly from the official projects,
review their current licenses and access conditions, and provide compatible local paths.
The official repositories are references; they do not guarantee that a checkpoint has
the exact project-specific classifier head expected by this pipeline.

用户必须直接从官方项目获取所需文件，阅读其当前许可证和访问条件，并提供兼容的
本地路径。官方原库链接仅用于指引；它们并不保证所下载权重包含本流程所需的特定分类头。

## Local path contract / 本地路径约定

Raw-WSI inference accepts command-line paths or these environment variables:

| Variable | Required content |
|---|---|
| `CRC_SURVIVAL_CELLVIT_ROOT` | local CellViT++ repository checkout |
| `CRC_SURVIVAL_CELLVIT_CHECKPOINT` | compatible CellViT++ segmentation checkpoint |
| `CRC_SURVIVAL_CELLVIT_CLASSIFIERS` | directory containing the seven compatible phenotype heads |
| `CRC_SURVIVAL_UNI_CHECKPOINT` | compatible UNI-derived coarse-classifier checkpoint |

All default locations are under the ignored `external/` directory. `vendor/`, `weights/`,
`external/`, and `*.pth` are ignored to prevent accidental publication.

默认路径均位于已忽略的 `external/` 目录下；`vendor/`、`weights/`、`external/`
及所有 `*.pth` 均被忽略，以防止误上传第三方资产。

## License boundary / 许可证边界

The subproject license applies only to repository-original code, documentation and the
included CRC survival-model files. It does not relicense any separately obtained UNI or
CellViT++ material. / 本子项目许可证仅覆盖仓库原创代码、文档及随附的 CRC 生存模型文件，
不对用户另行取得的 UNI 或 CellViT++ 材料重新授权。
