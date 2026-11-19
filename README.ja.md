<p align="center"><img src="web/public/favicon.svg" width="80" alt="mzDesk logo"></p>

<h1 align="center">mzDesk</h1>

<p align="center"><a href="README.md">English</a> | <a href="README.zh-CN.md">简体中文</a> | <b>日本語</b> | <a href="README.de.md">Deutsch</a> | <a href="README.fr.md">Français</a> | <a href="README.es.md">Español</a></p>

コード不要の低分子質量分析ツールです。中心となるワークフローは、MS/MS スペクトルから候補構造を順位付けし、それぞれにエビデンスレベルを付けます。

- **オンライン**：[mzdesk.com](https://mzdesk.com)（登録不要）
- **ローカル**：Windows、macOS、Linux でオフライン動作
- **ドキュメント**：[mzdesk.com/docs](https://mzdesk.com/docs)

![mzDesk によるカフェインのフラグメント解釈](.github/readme/screenshot.png)

## ローカル版のインストール

macOS と Linux：

```bash
curl -fsSL https://mzdesk.com/install.sh | sh
```

Windows（PowerShell）：

```powershell
irm https://mzdesk.com/install.ps1 | iex
```

管理者権限は不要です。インストール後、次のコマンドで起動します：

```bash
mzdesk
```

ブラウザで <http://localhost:8731> が開きます。計算エンジンとスペクトルライブラリは、ページで初めて必要になったときにダウンロードされます（合計約 9 GB）。[インストールの詳細](https://mzdesk.com/docs/local-version)

## 機能

![mzDesk：Web 版とローカル版、モジュールとワークフロー](.github/readme/overview.png)

## ライセンス

[MIT](LICENSE)。同梱およびダウンロードされるツールには、それぞれのライセンスが適用されます。[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) を参照してください。
