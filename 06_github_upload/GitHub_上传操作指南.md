# GitHub 仓库上传操作指南（geo-redundancy-sop）

> 本指南面向**不熟悉命令行的研究者**，提供三种上传方式，并把每一步都写清楚。
> 目标：把 `D:\AIwork\GEO\paper` 这个已整理好的投稿包，发布成一个公开的 GitHub 仓库
> `https://github.com/<你的用户名>/geo-redundancy-sop`，并在论文中引用它。

---

## 0. 当前已替你做好的事（重要）

你机器上的 `D:\AIwork\GEO\paper` 目录**已经是一个初始化好的本地 Git 仓库**，并且完成了第一次提交（commit），包含 30 个文件（稿件、报告、证据、清单、8 个脚本、6 张图）。所以你**不需要**再执行 `git init` / `git add` / `git commit`——只要创建远程仓库并推送即可（见方式一的「第 3–4 步」）。

> 如果你把文件复制到了别的位置，或想从零开始，请按方式一的「第 1–2 步」自行初始化。

---

## 1. 前置准备（只需做一次）

1. **注册 GitHub 账号**：https://github.com/signup （免费）。建议用户名用机构/姓名拼音，例如 `hebeu-fu`。
2. **安装 Git**：https://git-scm.com/downloads （Windows 版）。安装时一路默认即可，编辑器选 Notepad 或 VS Code 都行。
3. **配置本地身份**（打开 Git Bash 或终端，粘贴执行，仅一次）：
   ```bash
   git config --global user.name "Zexian Fu"
   git config --global user.email "fuzexian@hebeu.edu.cn"
   ```
4. **认证方式（二选一，推荐方式 A）**：
   - **A. Personal Access Token（HTTPS，最简单）**
     1. 登录 GitHub → 右上角头像 → **Settings → Developer settings → Personal access tokens → Tokens (classic) → Generate new token (classic)**。
     2. Note 填 `geo-redundancy-sop`；Expiration 选 `90 days` 或更长；勾选 `repo` 权限；点击 **Generate token**。
     3. **立刻复制那串 `ghp_xxx` 令牌并保存到安全地方**——页面关闭后不可再查看。
     4. 之后 `git push` 时，用户名填你的 GitHub 用户名，密码/口令**粘贴这个 token**（不是你的登录密码）。
   - **B. SSH 密钥（一劳永逸，略复杂）**
     1. 在 Git Bash 执行 `ssh-keygen -t ed25519 -C "fuzexian@hebeu.edu.cn"`，一路回车。
     2. 执行 `cat ~/.ssh/id_ed25519.pub`，复制输出内容。
     3. GitHub → Settings → **SSH and GPG keys → New SSH key**，粘贴保存。
     4. 之后远程地址用 `git@github.com:<用户名>/geo-redundancy-sop.git`（见第 3 步选 SSH 地址）。

---

## 2. 方式一：命令行（推荐，可复现、最稳妥）

### 第 1 步：在 GitHub 网页上创建「空仓库」
1. 登录 GitHub，点击右上角 **“+” → New repository**。
2. **Repository name** 填：`geo-redundancy-sop`
3. **Description** 填：`Reproducible SOP and tool for detecting cross-GSE sample redundancy and "face-swapped" datasets in colorectal cancer transcriptomics (GEO).`
4. 选 **Public**（公开，便于论文引用与审稿人核查）。
5. **不要**勾选 “Add a README file” / “Add .gitignore” / “Choose a license”——因为本地已经有了（避免冲突）。
6. 点击 **Create repository**。

### 第 2 步：复制新仓库的地址
创建后会进入一个快速上手页面。根据你选的认证方式复制地址：
- HTTPS：`https://github.com/<你的用户名>/geo-redundancy-sop.git`
- SSH：`git@github.com:<你的用户名>/geo-redundancy-sop.git`

### 第 3 步：关联远程仓库（在本地 paper 目录执行）
打开终端，`cd` 到 `D:\AIwork\GEO\paper`，然后执行（把地址换成你自己的）：
```bash
cd /d "D:\AIwork\GEO\paper"
git remote add origin https://github.com/<你的用户名>/geo-redundancy-sop.git
git branch -M main
```
> 如果提示 `remote origin already exists`，说明之前设过，可先 `git remote remove origin` 再重设，或 `git remote set-url origin <新地址>`。

### 第 4 步：推送
```bash
git push -u origin main
```
- 用 Token 认证时，弹出输入用户名/密码：用户名填 GitHub 账号名，**密码粘贴 token**。
- 成功后终端会显示进度条和 `main -> main`。

### 第 5 步：验证
浏览器打开 `https://github.com/<你的用户名>/geo-redundancy-sop`，应看到 30 个文件已在线。

---

## 3. 方式二：GitHub Desktop（图形界面，最省心）

适合完全不想碰命令行的用户。
1. 下载安装 https://desktop.github.com/ 并登录 GitHub 账号。
2. **File → Add Local Repository**，选择 `D:\AIwork\GEO\paper`（它已是 Git 仓库，能被识别）。
3. 点击 **Publish repository**（发布仓库）：
   - Name：`geo-redundancy-sop`
   - Description：同上
   - 勾选 **Keep this code private** 的**反选**（即公开发布）
   - 点 **Publish**
4. 之后若稿件有修改，在 GitHub Desktop 里写一句 Summary（如 `update manuscript with authors and funding`），点 **Commit to main**，再点 **Push origin** 即可同步。

---

## 4. 方式三：网页直接拖拽上传（仅应急，不推荐）

GitHub 仓库页面点 **Add file → Upload files**，把文件拖进去。
- **缺点**：每次最多上传有限文件、不能保留目录结构、大文件（如 `dedup_report.json` 2.9 MB）易失败、无法批量管理版本。**不建议用于正式发布。**

---

## 5. 推送后必做：把论文里的占位 URL 改成真实地址

稿件 `Manuscript_BIB_CRC_GEO_redundancy.md`（及导出的 .docx）中目前有两处占位：
```
https://github.com/[your-GitHub-username]/geo-redundancy-sop
```
请把它**全文替换**为你的真实地址，例如：
```
https://github.com/hebeu-fu/geo-redundancy-sop
```
改完后在本地 paper 目录重新提交并推送：
```bash
cd /d "D:\AIwork\GEO\paper"
# 用编辑器改好上面两处后：
git add -A
git commit -m "Update GitHub repository URL in manuscript"
git push
```
> 注意：`.docx` 是二进制，改 URL 需在 Word 里改完再覆盖文件，然后 `git add -A` 提交。

---

## 6. 让仓库更规范（加分项，约 10 分钟）

1. **加 LICENSE**：仓库页 → **Add file → Create new file**，文件名 `LICENSE`，GitHub 会自动提示选 MIT / Apache-2.0 模板（研究代码常用 **MIT** 或 **BSD-3-Clause**）。这能让别人合法复用。
2. **完善 README 顶部**：当前 `README.md` 已是中文清单，可在最上方加一句英文简介与引用方式，方便国际读者。
3. **加 Topics**：仓库页右侧 **About → gear 图标 → Topics**，添加 `colorectal-cancer`、`geo`、`meta-analysis`、`data-quality`、`bioinformatics`、`redundancy-detection` 等，提升可发现性。
4. **发 Release / 归档 DOI（可被正式引用）**：
   - 在仓库页 **Releases → Draft a new release**，填 `v1.0`，写一段说明，点 Publish。
   - 如需可引用 DOI，把仓库接入 **Zenodo**（https://zenodo.org，用 GitHub 账号登录 → GitHub → 选本仓库 → 开启 Archiving），每次打 Release 会自动生成一个 DOI（如 `10.5281/zenodo.XXXXXXX`）。论文中可写 “Code and data are available at GitHub (https://…) and archived at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX).”

---

## 7. 常见问题

- **推送被拒 `failed to push`**：通常是远程已有 README 与本地冲突。因为创建仓库时若勾了 “Add README”，需先 `git pull --rebase origin main` 再 `git push`。**本次指南要求创建时别勾选**，可避免此问题。
- **Token 被当密码拒绝**：确认粘贴的是 `ghp_` 开头的 token，不是登录密码；token 需有 `repo` 权限且未过期。
- **大文件警告**：本仓库最大文件约 2.9 MB，远未触及 GitHub 单文件 100 MB 限制，可正常推送。若日后加入 >50 MB 文件，建议用 Git LFS（本稿暂不需要）。
- **CRLF 换行警告**：Windows 下出现的 `LF will be replaced by CRLF` 是普通提示，不影响功能，可忽略；如需消除可在仓库根目录加 `.gitattributes` 文件，内容加一行 `* text=auto`。

---

## 8. 一句话流程总结

本地仓库已建好并提交 → 网页建同名空仓库（不勾 README）→ `git remote add origin <地址>` → `git branch -M main` → `git push -u origin main` → 改稿件占位 URL 再推一次 → 加 LICENSE/Topics/Release（可选）→ 论文引用真实地址。
