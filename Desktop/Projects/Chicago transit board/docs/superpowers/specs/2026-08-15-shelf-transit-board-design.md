# Shelf Transit Board — Design

**Date:** 2026-08-15
**Status:** Approved, ready for implementation planning
**Supersedes:** the CTA Live Wall Map build guide (`deepseek_latex_20260814_91d2ef.tex`)

---

## 1. What this is

A 10″ × 5″ amber LED board that sits on a bookshelf and shows how soon the next two
vehicles arrive on four specific trips:

| Row | Service | Origin | Direction |
|-----|---------|--------|-----------|
| 0 | CTA Red Line | Morse | southbound, toward 95th/Dan Ryan |
| 1 | CTA Purple Line | Howard | northbound, toward Linden |
| 2 | Metra UP-N | Rogers Park | inbound, toward Ogilvie (OTC) |
| 3 | CTA Bus 22 | Clark & Greenleaf | southbound |

It is WiFi-connected and has no controls. It shows numbers and nothing else.

### Relationship to the previous scope

The prior project was a wall-mounted printed CTA map with one addressable LED per
station — roughly 150 drilled holes, soldered strip segments, power injection, and a
station-to-LED-index CSV. That is abandoned. This shares only the choice of transit APIs.
The build is roughly a tenth of the work and about a twentieth of the wiring.

---

## 2. Architecture

```
   CTA Train Tracker ─┐
   CTA Bus Tracker   ─┼──► Cloudflare Worker ──► GET /board.json
   Metra GTFS-RT     ─┘    (holds all keys)          │
                                                     │ https, every 30 s
                                                     ▼
                                             ESP32 + MAX7219 chain
                                             (fetch and draw only)
```

The split exists for three reasons, in order of importance:

1. **Metra publishes GTFS-realtime protobuf covering every train on the system.** Decoding
   and filtering that on a microcontroller with ~300 KB of usable RAM is the single
   hardest part of the naive design. In a Worker it is a library call.
2. **API keys never live on the shelf device.**
3. **Route changes need no reflash.** Relabeling a row, swapping a stop, or changing
   direction is a Worker edit.

### 2.1 No cron, no KV

An earlier sketch had a 30-second cron trigger writing results to Workers KV. That does
not fit the free tier: KV allows on the order of 1,000 writes/day and a 30-second cron is
2,880. Instead the ESP32's own request drives the poll, and the Worker uses the Cache API
(`caches.default`), which carries no equivalent write cap.

Cache TTLs: **30 s** for the two CTA feeds, **60 s** for Metra — larger payload, and their
terms ask for gentler polling. At one device polling every 30 s this is ~2,880 requests/day
against CTA and ~1,440 against Metra, well inside the Workers free allowance of 100,000
requests/day.

---

## 3. The Worker

Single route, `GET /board.json`. Per request: read cache; on miss or stale, poll upstream,
filter to the four trips, write cache, return.

API keys are Worker secrets (`wrangler secret put`), never committed.

### 3.1 Response contract

Deliberately tiny and fixed-shape, so the firmware parser stays trivial:

```json
{
  "ts": 1755280000,
  "rows": [
    {"label": "RED", "mins": [4, 11], "state": "ok"},
    {"label": "PUR", "mins": [7, 22],  "state": "ok"},
    {"label": "UPN", "mins": [],       "state": "nosvc"},
    {"label": "22",  "mins": [3, 14],  "state": "ok"}
  ]
}
```

- `ts` — POSIX seconds, when the data was assembled (not when it was served from cache).
- `mins` — whole minutes until arrival, ascending, at most 2 entries. `0` means *due*.
- `state` — `ok` | `nodata` | `nosvc` | `error`, **per row**.

Rows are always present and always in the same order. A dead Metra feed must never blank
the Red Line row, so each upstream is fetched and failed independently.

### 3.2 Upstream endpoints

| Service | Endpoint | Format | Auth |
|---|---|---|---|
| CTA Train Tracker | `lapi.transitchicago.com/api/1.0/ttarrivals.aspx` | XML; `outputType=JSON` for JSON | key in query string |
| CTA Bus Tracker | `ctabustracker.com/bustime/api/v2/getpredictions` | `format=json` | separate key in query string |
| Metra | GTFS-realtime `TripUpdates` | protobuf | HTTP Basic |

The exact Metra hostname and path must be confirmed against their current developer
documentation during implementation rather than assumed from memory.

---

## 4. Firmware

ESP32 + `MD_Parola` / `MD_MAX72XX`, roughly 150 lines. Every 30 s: HTTPS GET, parse with
ArduinoJson, render four zones. It holds no API keys and no route logic — only the Worker
URL, WiFi credentials, and how to draw a row.

Zone mapping over 32 chained 8×8 matrices (8 four-in-one boards):

```cpp
P.setZone(0,  0,  7);   // row 0 — RED
P.setZone(1,  8, 15);   // row 1 — PUR
P.setZone(2, 16, 23);   // row 2 — UPN
P.setZone(3, 24, 31);   // row 3 — 22
```

Parola supports up to 16 zones on ESP32 and a zone is a run of consecutive modules, so this
needs no custom coordinate remapping.

`HARDWARE_TYPE` must match the boards. Most blue four-in-one boards are
`MD_MAX72XX::FC16_HW`, but a wrong value produces mirrored or scrambled text rather than a
blank display, which is confusing to diagnose. Confirm on first power-up.

---

## 5. Hardware

| # | Part | Qty | Note |
|---|---|---|---|
| 1 | MAX7219 four-in-one module, 32×8, FC-16 | 8 | all from one seller, one order |
| 2 | ESP32 DevKitC, 38-pin (WROOM-32) | 1 | plain ESP32, not S3, not ESP8266 |
| 3 | Female–female dupont jumpers | ~20 | 5 to the ESP32, 3 inter-row return cables |
| 4 | 5 V 2 A USB charger + cable | 1 | likely already owned |

Approximately $26 from AliExpress, $50 from Amazon. Cloudflare Worker: free tier, $0.

### 5.1 Wiring

Five wires from the ESP32 to the head of the chain, on the default VSPI pins:

```
VCC → 5V        DIN → GPIO23 (MOSI)
GND → GND       CLK → GPIO18 (SCK)
                CS  → GPIO5  (SS)
```

Boards are arranged two wide by four tall. Every row runs left to right, so three short
return cables route the chain from the right edge of one row back to the left edge of the
next. No serpentine layout — `MD_MAX72XX` does not remap it, and reversed rows would render
backwards.

### 5.2 Power

Measured draw is ~67 mA per 8×8 matrix at full brightness with a 62 kΩ Rset, scaling
linearly with the 16-step intensity register. 32 matrices is therefore ~2.1 A worst case,
but that assumes every LED lit at maximum. Four rows of sparse text at low intensity is
expected to be 200–400 mA.

Cheap clone boards frequently ship a lower Rset and draw well above the figure above, so:
**5 V 2 A charger, not 1 A**, and run a low intensity setting — which also looks better in
a room.

Feed 5 V to the module chain directly from the charger rather than through the ESP32's 5V
pin; 32 matrices is beyond what the dev board's USB path should carry. The ESP32 takes
power via its own USB connector. **Grounds must be common.**

---

## 6. Failure behaviour

A board confidently showing `4` when it does not know is worse than a board showing
nothing. Stale data must be visibly stale.

Every `state` value maps to exactly one rendering, so the firmware needs no judgement:

| `state` | Meaning | Row renders as |
|---|---|---|
| `ok` | predictions available | `LABEL  n  m` |
| `nodata` | feed responded, no predictions for this stop right now | `----` |
| `nosvc` | outside service hours for this trip | blank |
| `error` | upstream failed or timed out | `----` |

Conditions above the per-row level:

| Condition | Display |
|---|---|
| WiFi down, or Worker unreachable | all four rows `----` |
| Payload `ts` older than 3 minutes | all rows drop to lowest intensity |
| `mins` entry is `0` | that entry renders `DUE` |

Note the distinction between `----` (we do not know) and blank (we know nothing is
coming). `nodata` and `error` render identically on purpose — from the shelf, "the feed
has nothing for you" and "the feed is broken" are the same actionable fact, and
distinguishing them would cost characters the row does not have.

---

## 7. Correctness risks and how they are tested

Both of the following produce *plausible wrong numbers*, which is the dangerous failure
mode — a wrong answer that looks right invites no scrutiny.

### 7.1 Direction filtering

Morse serves the Red Line in both directions, Howard serves the Purple Line in both
directions, and Clark & Greenleaf has a stop on each side of the street. Selecting the
wrong one yields a perfectly reasonable two-digit number for a vehicle going the wrong way,
and nothing about the display would look wrong.

*Test:* recorded upstream responses as fixtures, asserting each row selects only the
intended direction. Fixtures must include a moment when both directions have imminent
arrivals, since that is when a direction bug is invisible in aggregate.

### 7.2 Timezone

CTA returns local wall-clock strings; Metra returns POSIX timestamps. The board is in
`America/Chicago` with DST. An hour error still renders a believable two-digit number.

*Test:* fixtures evaluated at a DST boundary, asserting computed minutes against
hand-checked expected values.

### 7.3 Stop ID resolution

Stop and route identifiers are **not** invented in this document. Morse and Howard `mapid`
values, the southbound Clark & Greenleaf `stpid`, and the Metra UP-N Rogers Park stop ID
are all resolved during implementation from GTFS static data and the Bus Tracker
`getstops` call, then recorded in the Worker config with a comment naming the source.

---

## 8. Risks

**Worker CPU limit on the Metra parse.** Free-tier Workers cap CPU time per request. A
systemwide GTFS-realtime TripUpdates parse may exceed it. The 60 s cache makes the parse
rare, which may be sufficient. If it is not, the fallback is to serve the Metra row from
**static GTFS schedule times** — UP-N Rogers Park inbound is roughly 40 departures per day,
small enough to bake into the Worker as a table. Scheduled rather than realtime is a
legitimate v1 for that one row; the other three are unaffected.

**API key lead time.** The CTA Bus Tracker key and Metra credentials are separate
applications that take days to approve. These are requested on day one, before any code is
written. The CTA Train Tracker key is quick by comparison.

---

## 9. Explicitly not doing

Night dimming schedule. Captive-portal WiFi setup — credentials are compile-time constants.
Enclosure or case. Web config UI. More than two upcoming arrivals per row. Service alerts
and delay reasons. Colour.

---

## 10. Deliverables

1. **Repo first** — Worker, firmware, tests, fixtures, README with BOM and wiring.

   ```
   chicago-transit-board/
     worker/
       board.js
       board.test.js
       fixtures/          recorded upstream responses
       wrangler.toml
     firmware/
       board.ino
       config.h           wifi credentials, worker URL
     README.md
   ```

2. **LaTeX build guide, generated after the physical build is verified working** — in the
   style of the previous wall-map document, with code listings inline. Written from a build
   known to run, not from the design.
