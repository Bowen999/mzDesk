import { useSyncExternalStore } from "react";

/**
 * Frontend-only password lock. Only the locked subgroups (DMPK, ICP-MS) sit behind it; every other module and every workflow is open.
 * This keeps casual viewers out; it is NOT access control (the API is open, and anyone with the
 * browser's dev tools can flip the flag), so never treat it as protection for real secrets.
 * Only a salted PBKDF2 hash is stored here; the password itself is not part of the source.
 */
const SALT_HEX = "7f4282fbfe80e0a447687b7fd2817b09";                                          // 16 random bytes
const HASH_HEX = "9bcbf15d44d26180b6d179732218c949367dae9a2e5be513211ce742c6d343b0";          // PBKDF2-HMAC-SHA256, 32 bytes
const ITERATIONS = 200_000;
const STORAGE_KEY = "mzdesk.unlocked";   // sessionStorage: cleared when the browser (tab) closes

/** Modules in these subgroups (DMPK, ICP-MS) are locked; every other module and every pipeline is open. */
export const LOCKED_SUBGROUPS = new Set(["dmpk", "inorganic"]);
export const isLockedModule = (m: { subgroup: string | null }) => !!m.subgroup && LOCKED_SUBGROUPS.has(m.subgroup);

const fromHex = (hex: string) => Uint8Array.from(hex.match(/../g) ?? [], b => parseInt(b, 16));

// --- unlock state: sessionStorage, with an in-memory fallback when storage is unavailable ---
let memory = false;
const listeners = new Set<() => void>();
const emit = () => listeners.forEach(l => l());

export function isUnlocked(): boolean {
  if (memory) return true;
  try { return sessionStorage.getItem(STORAGE_KEY) === "1"; } catch { return false; }
}

function setUnlocked(on: boolean) {
  memory = on;
  try { if (on) sessionStorage.setItem(STORAGE_KEY, "1"); else sessionStorage.removeItem(STORAGE_KEY); } catch { /* memory only */ }
  emit();
}

export function lock() { setUnlocked(false); }

/** Thrown when WebCrypto is missing (it only exists in secure contexts: https, localhost, 127.0.0.1). */
export class LockUnavailable extends Error {}

async function derive(password: string): Promise<Uint8Array> {
  const subtle = globalThis.crypto?.subtle;
  if (!subtle) throw new LockUnavailable("Password check needs a secure context. Open mzDesk via http://localhost or http://127.0.0.1.");
  const key = await subtle.importKey("raw", new TextEncoder().encode(password), "PBKDF2", false, ["deriveBits"]);
  const bits = await subtle.deriveBits({ name: "PBKDF2", hash: "SHA-256", salt: fromHex(SALT_HEX), iterations: ITERATIONS }, key, 256);
  return new Uint8Array(bits);
}

function sameBytes(a: Uint8Array, b: Uint8Array): boolean {
  let diff = a.length ^ b.length;
  for (let i = 0; i < Math.max(a.length, b.length); i++) diff |= (a[i] ?? 0) ^ (b[i] ?? 0);
  return diff === 0;
}

/** Checks the password; on success remembers it for this browser session. Rejects with LockUnavailable if WebCrypto is missing. */
export async function unlock(password: string): Promise<boolean> {
  const ok = sameBytes(await derive(password), fromHex(HASH_HEX));
  if (ok) setUnlocked(true);
  return ok;
}

export function subscribe(listener: () => void): () => void {
  listeners.add(listener);
  return () => { listeners.delete(listener); };
}

/** React hook: re-renders when the lock state changes. */
export const useUnlocked = () => useSyncExternalStore(subscribe, isUnlocked, () => false);
