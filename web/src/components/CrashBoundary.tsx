import { Component, type ErrorInfo, type ReactNode } from "react";
import { reportCrash } from "../api";
import { t } from "../i18n";
import { isLocalApp } from "../setup";

type Props = { children: ReactNode; label?: string };
type State = { error: Error | null; ref: string };

/**
 * Catches a page that fails to draw, so one broken page does not blank the whole app. The error, its stack and the page path
 * (no form data) are sent to the maintainer once (the local version keeps them on the computer: mzdesk/api/reports.py); the user gets the
 * reference and a way to try again.
 */
export class CrashBoundary extends Component<Props, State> {
  state: State = { error: null, ref: "" };

  static getDerivedStateFromError(error: Error): State {
    return { error, ref: "" };
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    void reportCrash(error, info.componentStack ?? "").then(ref => { if (ref) this.setState({ ref }); });
  }

  render() {
    const { error, ref } = this.state;
    if (!error) return this.props.children;
    // the reference code is set in <code>: it is a placeholder in the sentence, and the text around it is split off here
    const local = isLocalApp();
    const withRef = local ? t("The error was saved on this computer ({ref}). Your other pages are not affected.", { ref: "\u0001" })
      : t("The error was sent to the maintainer ({ref}). Your other pages are not affected.", { ref: "\u0001" });
    const [before, after = ""] = withRef.split("\u0001");
    return (
      <div className="page crash-page" role="alert">
        <h1>{this.props.label ? t("{label} stopped", { label: this.props.label }) : t("Something broke")}</h1>
        <p className="crash-text">{ref ? <>{before}<code>{ref}</code>{after}</> : (local ? t("The error was saved on this computer. Your other pages are not affected.")
          : t("The error was sent to the maintainer. Your other pages are not affected."))}</p>
        <div className="actions">
          <button type="button" className="btn" onClick={() => this.setState({ error: null, ref: "" })}>{t("Try again")}</button>
          <button type="button" className="btn secondary" onClick={() => window.location.reload()}>{t("Reload")}</button>
        </div>
      </div>
    );
  }
}
