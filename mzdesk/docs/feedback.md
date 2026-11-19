# Feedback

This page explains how to send feedback, what it contains, and what mzDesk records when something breaks.

## Send feedback

1. Click the bug button in the sidebar footer, next to the theme and the language.
2. Write what happened and what you were trying to do: 3 to 4000 characters.
3. Optional: add up to 3 images. Click **Add images**, drop files on the box, or paste a screenshot into the message.
4. Click **Send**.

The dialog then shows a reference such as `F-1A2B3C`. Quote it if you write again about the same thing.

In the local version, the bug button opens this form on mzdesk.com in a new tab. Say in your message that you use the local version, and on which page.

- Images can be PNG, JPG, WebP or GIF, each at most 20 MB. mzDesk scales them to at most 1600 pixels on the long side. Together they must stay under 5 MB after that.
- One address can send one message every 30 seconds, and 8 in an hour.

## What a message contains

| Sent | Not sent |
|---|---|
| Your message and images | Your inputs and results |
| The page you were on, such as `/modules/ms1_search` | The data of other pages |
| Your browser name and version | |
| The mzDesk version | |

The report does not store your IP address. If your inputs help, paste them into the message or add a screenshot.

## When something breaks

mzDesk files a crash report by itself in three cases. Each gets a reference such as `E-1A2B3C`.

| What you see | What the report holds |
|---|---|
| A page shows **… stopped**, with **Try again** and **Reload**. Your other pages keep working. | The error, its technical trace, the page and your browser. No inputs. |
| A run fails with a message that ends in `(ref E-…)`. | The error, its technical trace, and the names and sizes of the input fields. Never their values. |
| Nothing: a step of a workflow failed and the run went on. | The error and its technical trace. |

- The same error again within a day adds to the first report instead of making a new one.
- A message about your input, such as a missing field or an unreadable spectrum, is not a crash. Nothing is filed.

## Where it goes

### mzdesk.com

Messages and crash reports are stored on the mzDesk server, where the maintainer reads them. The server keeps the newest 500 reports.

### Local version

Your message goes through mzdesk.com, as above. Crash reports never leave your computer. They are kept in one file:

```text
~/.mzdesk/feedback/reports.sqlite
```

On Windows it is in `.mzdesk\feedback` in your user folder. Set `MZDESK_FEEDBACK_DIR` to use another folder. To tell the maintainer about a crash, send its reference with the bug button.

See [Privacy and Internet Search](/docs/privacy) for everything else that leaves your computer.
