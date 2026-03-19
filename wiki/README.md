# Wiki Content

This directory contains the markdown files that are synced to the [GitHub Wiki](https://github.com/hreinwald/drc/wiki).

## Files

- **Home.md** - The main landing page for the wiki
- **_Sidebar.md** - The sidebar navigation that appears on all wiki pages

## Syncing to GitHub Wiki

GitHub Wiki is a separate git repository. To sync the content from this directory to the actual GitHub Wiki:

### Option 1: Using the sync script (Recommended)

```bash
./sync-wiki.sh
```

This script will:
1. Clone the wiki repository to a temporary directory
2. Copy the wiki files from this directory
3. Commit and push the changes to the GitHub Wiki repository

### Option 2: Manual sync

```bash
# Clone the wiki repository
git clone https://github.com/hreinwald/drc.wiki.git

# Copy the files
cp wiki/Home.md drc.wiki/
cp wiki/_Sidebar.md drc.wiki/

# Commit and push
cd drc.wiki
git add .
git commit -m "Update wiki content"
git push origin master
```

## Important Notes

- The `wiki/` directory in the main repository is just a source directory
- Changes here won't automatically appear on GitHub Wiki
- You must sync the content using one of the methods above
- The GitHub Wiki is accessible at: https://github.com/hreinwald/drc/wiki
