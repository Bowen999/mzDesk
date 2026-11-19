<p align="center"><img src="web/public/favicon.svg" width="80" alt="mzDesk logo"></p>

<h1 align="center">mzDesk</h1>

<p align="center"><b>English</b> | <a href="README.zh-CN.md">简体中文</a> | <a href="README.ja.md">日本語</a> | <a href="README.de.md">Deutsch</a> | <a href="README.fr.md">Français</a> | <a href="README.es.md">Español</a></p>

Small-molecule mass spectrometry tools that need no code. The main workflow turns an MS/MS spectrum into ranked candidate structures, each with an evidence level.

- **Online:** [mzdesk.com](https://mzdesk.com), no sign-up
- **Local:** runs offline on Windows, macOS and Linux
- **Docs:** [mzdesk.com/docs](https://mzdesk.com/docs)

![Fragment explanation of caffeine in mzDesk](.github/readme/screenshot.png)

## Install the local version

macOS and Linux:

```bash
curl -fsSL https://mzdesk.com/install.sh | sh
```

Windows (PowerShell):

```powershell
irm https://mzdesk.com/install.ps1 | iex
```

No administrator rights needed. Then start it:

```bash
mzdesk
```

It opens <http://localhost:8731> in your browser. Engines and spectral libraries download the first time a page needs them (about 9 GB in total). [More on installing](https://mzdesk.com/docs/local-version)

## Functions

![mzDesk: the web and local versions, the modules and the workflows](.github/readme/overview.png)

## License

[MIT](LICENSE). Bundled and downloaded tools keep their own licenses: see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
