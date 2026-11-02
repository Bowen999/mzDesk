import assert from "node:assert/strict";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { after, test } from "node:test";
import { fileURLToPath } from "node:url";
import { createServer } from "vite";

const root = dirname(dirname(fileURLToPath(import.meta.url)));
const cache = mkdtempSync(join(tmpdir(), "mzdesk-network-test-"));
const server = await createServer({ configFile: false, root, cacheDir: cache, server: { middlewareMode: true, hmr: false }, appType: "custom" });
const { candidateNames, networkUses, isInternetOn, setInternet, setHostedSite } = await server.ssrLoadModule("/src/network.ts");

after(async () => { await server.close(); rmSync(cache, { recursive: true, force: true }); });

test("the public website has no switch: online lookups always run there, even after a 'No' saved by the local version", () => {
  setInternet(false);
  assert.equal(isInternetOn(), false);
  setHostedSite(true);
  assert.equal(isInternetOn(), true);
  setHostedSite(false);
  assert.equal(isInternetOn(), false, "the local version keeps the person's choice");
  setInternet(true);
});

test("a candidate list asks for the internet only for lines that are names", () => {
  assert.deepEqual(candidateNames("CCO\nc1ccccc1\nC(C(=O)O)N\n[NH3+]CC(=O)[O-]\nCl[Na]\nInChI=1S/CH4/h1H4"), []);
  assert.deepEqual(candidateNames("alanine\n\n# a comment\n  L-lactic acid \nCCO\nUrea\n56-12-2\nCo"), ["alanine", "L-lactic acid", "Urea", "56-12-2", "Co"]);
});

test("a line that only looks like a SMILES is a name when its rings or brackets cannot close", () => {
  assert.deepEqual(candidateNames("F6P\nC1\nPCB-187\nCI-1018\nPC(16:0/18:1)\nC%10\nCC(C\nCC)C"), ["F6P", "C1", "PCB-187", "CI-1018", "PC(16:0/18:1)", "C%10", "CC(C", "CC)C"]);
  assert.deepEqual(candidateNames("C1CC1\nc1ccc2ccccc2c1\nC%10CCCCC%10\n[13CH4]\n[NH3+]CC(=O)[O-]\nC1CC1C1CC1\nO=C1CCC(=O)N1"), []);       // digits inside [..] are not ring closures
});

test("a line with a space is a name, as the backend reads it", () => {
  assert.deepEqual(candidateNames("PC 33:1\nCC(=O)O acetic acid\nCCO"), ["PC 33:1", "CC(=O)O acetic acid"]);
});

const field = { type: "textarea", id: "candidates", label: "Candidate structures", network_names: true };
test("networkUses names the lookup when the candidate field holds a name", () => {
  assert.deepEqual(networkUses([field], { candidates: "CCO\nC(C(=O)O)N" }), []);
  assert.deepEqual(networkUses([field], { candidates: "" }), []);
  assert.deepEqual(networkUses([field], { candidates: "alanine" }), ["A candidate name is looked up in PubChem and UniChem."]);
  assert.deepEqual(networkUses([field], { candidates: "alanine\nserine\nCCO" }), ["Candidate names are looked up in PubChem and UniChem."]);
  assert.deepEqual(networkUses([{ ...field, network_names: undefined }], { candidates: "alanine" }), []);       // a textarea without the flag never asks
});
