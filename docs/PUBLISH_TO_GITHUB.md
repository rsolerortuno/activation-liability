# Publishing the complete repository to GitHub

The repository must be populated from the complete verified release, not by uploading only the
README and top-level documentation.

## Required GitHub permissions

The account or token performing the push needs:

- repository **Contents: read and write**;
- **Workflows/Actions: read and write**, because `.github/workflows/ci.yml` is versioned;
- access to `rsolerortuno/activation-liability`.

## Verified source archive

- Drive file ID: `1obHMFbJhyQ_9JDRstEt2ogpkmGm84uBa`
- Expected SHA-256: recorded in the adjacent `.zip.sha256` file.

Use `notebooks/publish_complete_repository_to_github.ipynb` for a guided Colab publication or
replace the current worktree locally with the extracted archive and push it with Git.
