<p align="center"><img src="web/public/favicon.svg" width="80" alt="mzDesk logo"></p>

<h1 align="center">mzDesk</h1>

<p align="center"><a href="README.md">English</a> | <a href="README.zh-CN.md">简体中文</a> | <a href="README.ja.md">日本語</a> | <a href="README.de.md">Deutsch</a> | <a href="README.fr.md">Français</a> | <b>Español</b></p>

Herramientas de espectrometría de masas de moléculas pequeñas, sin programar. El flujo de trabajo principal convierte un espectro MS/MS en una lista ordenada de estructuras candidatas, cada una con un nivel de evidencia.

- **En línea:** [mzdesk.com](https://mzdesk.com), sin registro
- **Local:** funciona sin conexión en Windows, macOS y Linux
- **Documentación:** [mzdesk.com/docs](https://mzdesk.com/docs)

![Explicación de fragmentos de la cafeína en mzDesk](.github/readme/screenshot.png)

## Instalar la versión local

macOS y Linux:

```bash
curl -fsSL https://mzdesk.com/install.sh | sh
```

Windows (PowerShell):

```powershell
irm https://mzdesk.com/install.ps1 | iex
```

No requiere permisos de administrador. Después, inícielo:

```bash
mzdesk
```

Se abre <http://localhost:8731> en el navegador. Los motores de cálculo y las bibliotecas de espectros se descargan la primera vez que una página los necesita (unos 9 GB en total). [Más sobre la instalación](https://mzdesk.com/docs/local-version)

## Funciones

![mzDesk: versiones web y local, módulos y flujos de trabajo](.github/readme/overview.png)

## Licencia

[MIT](LICENSE). Las herramientas incluidas o descargadas conservan sus propias licencias: consulte [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
