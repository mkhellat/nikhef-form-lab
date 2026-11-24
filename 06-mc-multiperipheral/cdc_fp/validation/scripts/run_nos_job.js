#!/usr/bin/env node
/*
 * Submit validation/ftn/cdc_fp_cases.f to a NOS 2.8.7 DtCyber and
 * write the printer listing.
 *
 * Run with cwd irrelevant. Requires:
 *   DTCYBER_NOS  NOS2.8.7 directory (cyber.ini, disks/)
 *   DTCYBER_ROOT DtCyber source root that contains dtcyber and automation/
 * Optional:
 *   CDC_JOB_DECK path of the generated deck (default: decks/cdcfpc.job
 *                under DTCYBER_NOS, not deleted)
 *
 * The Fortran file is argv[2] if given, else ../ftn/cdc_fp_cases.f
 * relative to this script. The listing is argv[3] if given, else stdout.
 *
 * Deadstart uses the [cyber] profile. Telnet is not required. The card
 * reader is channel 12 equipment 4, as in NOS2.8.7/cyber.ini.
 *
 * The dayfile on NOS 2.8.7 871 ended at STOP and did not print
 * "*** CDCFPC COMPLETE". After three minutes the script falls back to
 * the printer file LP5xx_C12_E5 when that file contains CDCVAL.
 */

const fs = require("fs");
const path = require("path");

const nosDir = process.env.DTCYBER_NOS;
const root = process.env.DTCYBER_ROOT;
if (!nosDir || !root) {
  console.error("DTCYBER_NOS and DTCYBER_ROOT are required");
  process.exit(1);
}

const DtCyber = require(path.join(root, "automation", "DtCyber"));
const ftnPath = process.argv[2]
  || path.join(__dirname, "..", "ftn", "cdc_fp_cases.f");
const listingPath = process.argv[3] || "";

const source = fs.readFileSync(ftnPath, "utf8").replace(/\r\n/g, "\n");
const deck = [
  "CDCFPC.",
  "USER,GUEST,GUEST.",
  "FTN5,L=0,OPT=0,GO.",
  "~eor",
  source.replace(/\n$/, ""),
  "~eor",
  "$PACK,OUTPUT.",
  "*** CDCFPC COMPLETE",
  "~eoi",
  "",
].join("\n");

process.chdir(nosDir);
fs.mkdirSync("decks", { recursive: true });
const deckPath = path.join(nosDir, "decks", "cdcfpc.job");
fs.writeFileSync(deckPath, deck);

const dtc = new DtCyber();
let listing = "";

function withTimeout(promise, ms) {
  let timer;
  const timeout = new Promise((_, reject) => {
    timer = setTimeout(() => reject(new Error("job wait timed out")), ms);
  });
  return Promise.race([promise, timeout]).finally(() => clearTimeout(timer));
}

function capturedListing() {
  if (listing.includes("CDCVAL")) return listing;
  const printer = path.join(nosDir, "LP5xx_C12_E5");
  try {
    const disk = fs.readFileSync(printer, "latin1");
    if (disk.includes("CDCVAL")) return disk;
  } catch (e) {
    /* printer file is not there yet */
  }
  return listing;
}

function writeListing(text) {
  if (listingPath) fs.writeFileSync(listingPath, text);
  else process.stdout.write(text);
}

dtc.start({
  detached: true,
  /* Keep stdin open. DtCyber treats EOF on stdin as shutdown, and a
     closed parent stdin (or /dev/null) delivers that EOF immediately. */
  stdio: ["pipe", "ignore", 2],
  unref: false,
})
.then(() => dtc.sleep(5000))
.then(() => dtc.connect())
.then(() => dtc.expect([{ re: /Operator> $/ }]))
.then(() => dtc.attachPrinter("LP5xx_C12_E5"))
.then(() => dtc.expect([{ re: /QUEUE FILE UTILITY COMPLETE/ }], "printer"))
.then(() => dtc.say("Deadstart complete, submitting CDCFPC"))
.then(() => withTimeout(dtc.runJob(12, 4, deckPath, data => {
  listing += data;
}), 180000))
.then(() => {
  const text = capturedListing();
  if (!text.includes("CDCVAL")) {
    throw new Error("job ended without CDCVAL lines");
  }
  writeListing(text);
  return dtc.shutdown(false);
})
.then(() => {
  process.exit(0);
})
.catch(err => {
  const text = capturedListing();
  const timedOut = err && String(err.message || err).includes("timed out");
  if (timedOut && text.includes("CDCVAL")) {
    writeListing(text);
    console.error(
      "job marker missing; wrote printer listing "
      + (listingPath || "stdout")
    );
    dtc.shutdown(false).catch(() => {}).then(() => process.exit(0));
    return;
  }
  if (listingPath && text) {
    try { fs.writeFileSync(listingPath, text); } catch (e) { /* keep the original error */ }
  } else if (text) {
    process.stderr.write(text);
  }
  console.error(err);
  dtc.shutdown(false).catch(() => {}).then(() => process.exit(1));
});
