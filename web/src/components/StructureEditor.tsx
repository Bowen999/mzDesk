import { useRef } from "react";
import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import "ketcher-react/dist/index.css";

const structServiceProvider = new StandaloneStructServiceProvider();

type Props = {
  smiles: string;                                   // 进入 Draw 时画布的初始结构
  molfile?: string;                                 // 切走再切回时恢复画布
  onChange: (v: { value: string; molfile: string }) => void;
};

/** Ketcher 画布（浏览器内运行，不联网）。每次改动后把 SMILES 和 molfile 交回表单。 */
export default function StructureEditor({ smiles, molfile, onChange }: Props) {
  const ketcher = useRef<any>(null);

  const sync = async () => {
    const k = ketcher.current;
    if (!k) return;
    const mol = await k.getMolfile();
    // 空画布的 molfile 只有 0 个原子
    const empty = /^\s*0\s+0\s/m.test(mol.split("\n")[3] ?? "");
    onChange(empty ? { value: "", molfile: "" } : { value: await k.getSmiles(), molfile: mol });
  };

  return (
    <div className="structure-editor">
      <Editor
        staticResourcesUrl=""
        structServiceProvider={structServiceProvider as any}
        errorHandler={(m: string) => console.warn("Ketcher:", m)}
        onInit={async (k: any) => {
          ketcher.current = k;
          const start = molfile || smiles;
          if (start) { try { await k.setMolecule(start); } catch { /* 起始结构读不进就留空画布 */ } }
          k.editor.subscribe("change", () => { void sync(); });
        }}
      />
    </div>
  );
}
