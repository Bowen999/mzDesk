import { Component, type ErrorInfo, type ReactNode } from "react";
import { reportCrash } from "../api";

type Props = { children: ReactNode; label?: string };
type State = { error: Error | null; ref: string };

/**
 * Catches a page that fails to draw, so one broken page does not blank the whole app. The error, its stack and the page path
 * (no form data) are sent to the maintainer once; the user gets the reference and a way to try again.
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
    return (
      <div className="page crash-page" role="alert">
        <h1>{this.props.label ? `${this.props.label} stopped` : "Something broke"}</h1>
        <p className="crash-text">The error was sent to the maintainer{ref ? <> (<code>{ref}</code>)</> : ""}. Your other pages are not affected.</p>
        <div className="actions">
          <button type="button" className="btn" onClick={() => this.setState({ error: null, ref: "" })}>Try again</button>
          <button type="button" className="btn secondary" onClick={() => window.location.reload()}>Reload</button>
        </div>
      </div>
    );
  }
}
