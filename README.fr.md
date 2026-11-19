<p align="center"><img src="web/public/favicon.svg" width="80" alt="mzDesk logo"></p>

<h1 align="center">mzDesk</h1>

<p align="center"><a href="README.md">English</a> | <a href="README.zh-CN.md">简体中文</a> | <a href="README.ja.md">日本語</a> | <a href="README.de.md">Deutsch</a> | <b>Français</b> | <a href="README.es.md">Español</a></p>

Des outils de spectrométrie de masse des petites molécules, sans programmation. Le workflow principal transforme un spectre MS/MS en une liste classée de structures candidates, chacune avec un niveau de preuve.

- **En ligne :** [mzdesk.com](https://mzdesk.com), sans inscription
- **En local :** fonctionne hors ligne sous Windows, macOS et Linux
- **Documentation :** [mzdesk.com/docs](https://mzdesk.com/docs)

![Explication des fragments de la caféine dans mzDesk](.github/readme/screenshot.png)

## Installer la version locale

macOS et Linux :

```bash
curl -fsSL https://mzdesk.com/install.sh | sh
```

Windows (PowerShell) :

```powershell
irm https://mzdesk.com/install.ps1 | iex
```

Aucun droit d'administrateur n'est nécessaire. Puis lancez :

```bash
mzdesk
```

La page <http://localhost:8731> s'ouvre dans votre navigateur. Les moteurs de calcul et les bibliothèques de spectres sont téléchargés la première fois qu'une page en a besoin (environ 9 Go au total). [En savoir plus sur l'installation](https://mzdesk.com/docs/local-version)

## Fonctions

![mzDesk : versions web et locale, modules et workflows](.github/readme/overview.png)

## Licence

[MIT](LICENSE). Les outils fournis ou téléchargés conservent leurs propres licences : voir [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
