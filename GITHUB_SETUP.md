# GitHub Setup Guide

Follow these steps to create a GitHub repository for your ArXiv summary tool.

## Step 1: Create Repository on GitHub

1. Go to https://github.com/new
2. Fill in the details:
   - **Repository name**: `arxiv-soft-matter-summary`
   - **Description**: "Daily summaries of arXiv soft condensed matter papers"
   - **Public** or **Private** (your choice)
   - **DO NOT** check "Initialize with README" (we already have one)
3. Click "Create repository"

## Step 2: Initialize Git Locally

Open PowerShell or Command Prompt in your project folder and run:

```bash
git init
git add .
git commit -m "Initial commit: ArXiv soft matter daily summary tool"
```

## Step 3: Connect to GitHub

GitHub will show you commands after creating the repo. Use these:

```bash
git remote add origin https://github.com/YOUR_USERNAME/arxiv-soft-matter-summary.git
git branch -M main
git push -u origin main
```

Replace `YOUR_USERNAME` with your actual GitHub username.

## Step 4: Verify

Go to your repository URL and you should see all your files!

## Alternative: Using GitHub Desktop

1. Download GitHub Desktop: https://desktop.github.com/
2. Open GitHub Desktop
3. File → Add Local Repository → Select your folder
4. Click "Publish repository"
5. Choose name and public/private setting
6. Click "Publish"

## Future Updates

After making changes to your code:

```bash
git add .
git commit -m "Description of changes"
git push
```

## Adding a Screenshot (Optional but Recommended)

1. Run the script once to generate `arxiv_summary.html`
2. Take a screenshot of the HTML output
3. Save as `screenshot.png` in your project folder
4. Update git:
   ```bash
   git add screenshot.png
   git commit -m "Add example screenshot"
   git push
   ```

## Common Issues

**"Git not recognized"**
- Install Git: https://git-scm.com/downloads
- Restart PowerShell after installation

**Authentication required**
- Use GitHub Personal Access Token instead of password
- Settings → Developer settings → Personal access tokens → Generate new token
- Use token as password when pushing

**Or use SSH:**
```bash
git remote set-url origin git@github.com:YOUR_USERNAME/arxiv-soft-matter-summary.git
```

## That's it!

Your code is now on GitHub and ready to share! 🎉
