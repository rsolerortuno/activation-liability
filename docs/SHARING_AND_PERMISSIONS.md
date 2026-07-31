# Sharing and permissions

## Recommended public configuration

For the Google Drive root (https://drive.google.com/drive/folders/1f3K-MzEQDsUFmMIb5cnFxH47_G-HYPKN):

1. Set **General access** to **Anyone with the link**.
2. Select **Viewer**, not Editor.
3. Keep the owner account and named collaborators as Editors.
4. Do not expose credentials, API keys, controlled-access data or personal information.

## Can the assistant continue editing?

Yes, provided the Google account connected to this conversation still has edit permission. Public
link settings do not grant or remove the connected account's owner/editor rights. Making everyone
with the link an Editor is unnecessary and unsafe: anyone could overwrite, move or delete files.

## Integrity controls

- Publish SHA-256 files for every ZIP release.
- Keep source and result manifests under version control.
- Treat Drive as a data/release mirror, not the authoritative source-code history.
- Use GitHub for reviewed code changes and immutable release tags.
