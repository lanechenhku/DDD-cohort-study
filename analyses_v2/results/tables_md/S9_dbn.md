**Table S9a. Feature-level DBN: arcs with bootstrap strength ≥ 0.3.**

| From (scan *i*) | To (scan *i*+1) | Type | In full-data network | Bootstrap strength |
|---|---|---|---|---|
| SS | DB | cross | yes | 1.000 |
| DN | MC | cross | yes | 1.000 |
| MC | DN | cross | yes | 1.000 |
| DB | DN | cross | yes | 1.000 |
| DB | DH | cross | yes | 0.935 |
| DB | MC | cross | yes | 0.905 |
| SS | HIZ | cross | yes | 0.790 |
| DB | HIZ | cross | yes | 0.740 |
| DB | SS | cross | yes | 0.500 |
| MC | MC | autoregressive | yes | 1.000 |
| DB | DB | autoregressive | yes | 1.000 |
| SN | SN | autoregressive | yes | 1.000 |
| SS | SS | autoregressive | yes | 1.000 |
| DN | DN | autoregressive | yes | 1.000 |
| HIZ | HIZ | autoregressive | yes | 1.000 |
| DH | DH | autoregressive | yes | 0.715 |

**Table S9b. Disc-level DBN: arcs with bootstrap strength ≥ 0.3.**

| From (scan *i*) | To (scan *i*+1) | Type | In full-data network | Bootstrap strength |
|---|---|---|---|---|
| L3 | L2 | cross | yes | 1.000 |
| L4 | L3 | cross | yes | 1.000 |
| L3 | L4 | cross | yes | 0.860 |
| L3 | L1 | cross | yes | 0.720 |
| L1 | L2 | cross | no | 0.540 |
| L2 | L1 | cross | no | 0.485 |
| L1 | L3 | cross | no | 0.315 |
| L4 | L4 | autoregressive | yes | 1.000 |
| L1 | L1 | autoregressive | yes | 1.000 |
| L5 | L5 | autoregressive | yes | 1.000 |
| L3 | L3 | autoregressive | yes | 1.000 |
| L2 | L2 | autoregressive | yes | 1.000 |
