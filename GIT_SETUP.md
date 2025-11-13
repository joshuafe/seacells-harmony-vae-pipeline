# Git Setup and GitHub Push Guide

This guide will help you initialize the repository and push to GitHub.

## Step 1: Initialize Git Repository

```bash
# Initialize git in the current directory
git init

# Check git status (see what will be committed)
git status
```

## Step 2: Add Files to Git

```bash
# Add all files (respects .gitignore)
git add .

# Review what will be committed
git status
```

You should see files added (green) and ignored files NOT shown (data/, outputs/, etc.)

## Step 3: Create Initial Commit

```bash
# Create the first commit
git commit -m "Initial commit: SEACells + Harmony + VAE pipeline

- Per-sample SEACells with configurable compression ratio
- Harmony batch correction for sample integration
- PhenoGraph and Leiden clustering analysis
- VAE training with Optuna optimization
- Docker support for reproducibility
- Comprehensive visualization scripts
"
```

## Step 4: Create GitHub Repository

### Option A: Using GitHub CLI (gh)

```bash
# Login to GitHub (if not already)
gh auth login

# Create a new repository
gh repo create seacells-harmony-vae-pipeline --public --source=. --remote=origin

# Push to GitHub
git push -u origin main
```

### Option B: Using GitHub Web Interface

1. Go to https://github.com/new
2. Repository name: `seacells-harmony-vae-pipeline`
3. Description: "Flow cytometry analysis pipeline with SEACells, Harmony batch correction, clustering, and VAE"
4. Choose Public or Private
5. **DO NOT** initialize with README, .gitignore, or license (we already have these)
6. Click "Create repository"

Then connect and push:

```bash
# Add the remote (replace YOUR_USERNAME)
git remote add origin https://github.com/YOUR_USERNAME/seacells-harmony-vae-pipeline.git

# Verify remote
git remote -v

# Push to GitHub
git branch -M main
git push -u origin main
```

## Step 5: Verify on GitHub

Visit your repository on GitHub and verify:
- ✓ All code files are present
- ✓ Data files are NOT present (ignored)
- ✓ Output directories are NOT present (ignored)
- ✓ README.md displays properly

## Common Git Workflow Going Forward

### Make Changes and Commit

```bash
# See what changed
git status
git diff

# Add specific files
git add scripts/my_script.py

# Or add all changes
git add .

# Commit with descriptive message
git commit -m "Add feature X to do Y"

# Push to GitHub
git push
```

### Before Making Changes (Pull Latest)

```bash
# Always pull before starting work
git pull
```

### Create Feature Branches

```bash
# Create and switch to new branch
git checkout -b feature/new-analysis

# Make changes and commit
git add .
git commit -m "Add new analysis feature"

# Push branch to GitHub
git push -u origin feature/new-analysis

# Create pull request on GitHub
gh pr create --title "Add new analysis feature" --body "Description of changes"
```

### Test Different Compression Ratios (Using Branches)

```bash
# Create branch for 30:1 compression
git checkout -b experiment/compression-30

# Modify or just run with different settings
python run_pipeline.py --compression-ratio 30.0

# Commit results (if you want to track)
git add .
git commit -m "Test 30:1 compression ratio"
git push -u origin experiment/compression-30

# Switch back to main
git checkout main

# Try another ratio
git checkout -b experiment/compression-75
python run_pipeline.py --compression-ratio 75.0
```

## Useful Git Commands

```bash
# View commit history
git log --oneline

# View changes since last commit
git diff

# Undo changes to a file (before commit)
git checkout -- filename

# View all branches
git branch -a

# Switch branches
git checkout branch-name

# Delete local branch
git branch -d branch-name

# View remote info
git remote -v

# See what's ignored
git status --ignored
```

## Important Notes

### What IS Tracked (In Git):
- ✓ All Python scripts
- ✓ Configuration files (environment.yml, Dockerfile)
- ✓ Pipeline orchestration (run_pipeline.py, run_pipeline.sh)
- ✓ Documentation (README.md, guides)
- ✓ SEACells package code

### What is NOT Tracked (Ignored):
- ✗ Data files (*.fcs, *.csv, *.h5ad)
- ✗ Output directories (seacells_output/, metacell_vae_output/)
- ✗ Plots and visualizations (*.png, *.pdf)
- ✗ Model files (*.pt, *.db)
- ✗ Log files (*.log)
- ✗ Python cache (__pycache__, *.pyc)

This keeps the repository clean and focused on code, not data or results.

## Collaborating with Others

### Clone Repository

```bash
# Clone to a new location
git clone https://github.com/YOUR_USERNAME/seacells-harmony-vae-pipeline.git
cd seacells-harmony-vae-pipeline

# Set up environment
conda env create -f environment.yml
conda activate metacell_vae_env
cd SEACells && pip install -e . && cd ..
```

### Pull Request Workflow

1. Fork the repository (if not owner)
2. Create feature branch: `git checkout -b feature/my-feature`
3. Make changes and commit
4. Push to GitHub: `git push -u origin feature/my-feature`
5. Create pull request on GitHub
6. Review and merge

## Tagging Releases

When you have a stable version:

```bash
# Create annotated tag
git tag -a v1.0.0 -m "First stable release with configurable compression"

# Push tags to GitHub
git push --tags

# View tags
git tag -l
```

## Configuration

### Set Your Identity (First Time Only)

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### View Configuration

```bash
git config --list
```

## Getting Help

```bash
# Help for any command
git help <command>
git help commit
git help push
```

## Ready to Commit?

Run these commands now:

```bash
# 1. Initialize
git init

# 2. Add files
git add .

# 3. Check what will be committed
git status

# 4. Create initial commit
git commit -m "Initial commit: SEACells + Harmony + VAE pipeline"

# 5. Create GitHub repo and push (see Option A or B above)
```

---

**After pushing to GitHub, you can safely experiment with different compression ratios knowing your code is backed up!**
