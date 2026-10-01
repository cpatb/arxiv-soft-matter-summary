"""Silver to Color studio: turn a restored film plus a script into finished videos.

Stages, each a module:

    colorize   B&W source -> decimated, deflickered B&W master + aligned colour master
               (OpenVINO on the Intel iGPU, ONNX Runtime CPU fallback)
    episode    YAML episode spec -> validated cue list
    voice      cue text -> narration audio + word timings (edge-tts scratch voice,
               or your own recordings, which always win)
    timeline   cues + audio durations -> absolute timeline, chapters, captions
    render     timeline -> full 16:9 video, vertical Shorts, SRT, thumbnails

Run `python -m studio --help`.
"""

import sys
from pathlib import Path

# The restoration helpers (temporal_chroma, pipeline, make_metadata) live in ../tools.
_TOOLS = Path(__file__).resolve().parent.parent / "tools"
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

__all__ = ["colorize", "episode", "voice", "timeline", "render"]
