# Session Report

This page explains the Session report: one download with the results of every page you ran in this browser tab.

## Download it

1. Run one or more pages.
2. Click **Session report** in the sidebar footer. In the folded sidebar it is an icon.
3. Your browser saves a file such as `mzdesk-20261031-153012-MDT-k7q2m.zip`.

- The button is disabled until a page has a result.
- The number on the button is how many pages the report holds.
- The name is the run ID: the date and time, the time zone and five random characters.

## Which pages it holds

Every module or workflow page that has a result in this tab, in the order they finished.

A page leaves the report when you:

- click **Clear**, or load an example;
- run it again (it comes back when the new run finishes);
- run it and the run fails.

Moving to another page keeps its result. Reloading the tab or closing it empties the report.

## What is in the file

The `.zip` holds one folder with two files.

| File | For | Contains |
|---|---|---|
| `report.html` | Reading and printing | Every page as it looks on screen, in one file that opens offline. |
| `report.json` | Re-use and checking | The inputs and the raw results of every page, and the versions of the tools. |

### report.html

- The top lists the run ID, the time, the mzDesk version and the pages.
- Each page starts with its inputs, as short chips, and the inputs as sent (JSON, folded).
- Then the same results as on screen: the answer, tables, figures and, for a workflow, the steps.
- Folded blocks open with a click. Buttons that need the app, such as filters, table pages and **View**, are left out.
- A table shows at most 500 rows. Long input text is cut after 1500 characters. `report.json` has everything.
- Figures are pictures: no zoom or hover.
- Each page starts on a new sheet when you print.
- The report uses the theme you have on screen, light or dark.

### report.json

| Key | Meaning |
|---|---|
| `run_id`, `generated_at` | The run ID and the time the report was made. |
| `app` | mzDesk and its version. |
| `environment` | The versions of Python, the packages, the engines, the models and the data. |
| `pages` | One entry per page. |

Each entry of `pages` has the page `id` and `title`, `finished_at`, `seconds`, `inputs`, `result`, `steps` (workflows), and `based_on`: the tools the page relies on, each with its version.

## Versions

The versions come from the mzDesk that ran the pages: the server of mzdesk.com, or your own computer for the [local version](/docs/local-version).

- A version of `null` means the tool is not installed there, or its version could not be found.
- An online service, such as PubChem, has no version. Its entry records when it was asked.
- If the versions cannot be read within 20 seconds, `environment` is `null`, `environment_error` says why, and the report is made anyway.

## Times

All times in the report are Edmonton time (`America/Edmonton`), whatever your computer's time zone. `report.json` gives the zone and the UTC offset.

## Where it is made

Your browser builds the report from the results on screen and saves it. Nothing is uploaded to make it. The only request is for the versions of the tools. See [Privacy and Internet Search](/docs/privacy).

If one page cannot be drawn, `report.html` says so at its place, and its result is still in `report.json`.
