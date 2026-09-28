"""Telling a person something happened, from a TUI they may not be looking at.

Over ssh the only channel back to the person's desktop is the terminal
stream itself, so an event goes out as escape sequences their terminal may
turn into a notification: the bell, OSC 9 (iTerm2, WezTerm, Windows
Terminal, ...) and OSC 777 (urxvt, foot, Ghostty, ...). A terminal ignores
the ones it does not know. Inside tmux each is wrapped in tmux's passthrough
envelope, which tmux only honours with `set -g allow-passthrough on`.

For anything further -- a phone, an email -- `--on-event CMD` runs a command
per event with the event in its environment (ntfy, mail, a webhook)."""
from __future__ import annotations

import os
import subprocess

from ..core.events import Event

NOTIFY_MODES = ("off", "bell", "osc", "all")


def tmux_wrap(sequence: str) -> str:
    """tmux passthrough: DCS `tmux;`, the sequence with every ESC doubled,
    then ST."""
    return "\x1bPtmux;" + sequence.replace("\x1b", "\x1b\x1b") + "\x1b\\"


def _clean(text: str) -> str:
    # Neither a BEL nor an ESC may appear inside an OSC payload: either would
    # end it early and dump the rest onto the screen.
    return "".join(c for c in text if c >= " " and c != "\x7f")


def sequences(event: Event, mode: str, in_tmux: bool | None = None) -> str:
    """What to write to the terminal for `event` under `mode`."""
    if mode == "off":
        return ""
    in_tmux = bool(os.environ.get("TMUX")) if in_tmux is None else in_tmux
    parts = []
    if mode in ("bell", "all"):
        parts.append("\a")
    if mode in ("osc", "all"):
        text = _clean(event.text)
        for seq in (f"\x1b]9;{text}\x07", f"\x1b]777;notify;orcamon;{text}\x07"):
            parts.append(tmux_wrap(seq) if in_tmux else seq)
    return "".join(parts)


class Hook:
    """`--on-event CMD`: run CMD through the shell once per event, without
    waiting for it, with its output discarded. Children are reaped on later
    calls so none lingers as a zombie."""

    def __init__(self, command: str | None):
        self.command = command
        self._children: list[subprocess.Popen] = []

    def __call__(self, event: Event) -> None:
        self._children = [p for p in self._children if p.poll() is None]
        if not self.command:
            return
        env = {
            **os.environ,
            "ORCAMON_EVENT": event.kind,
            "ORCAMON_JOB": event.job,
            "ORCAMON_STATUS": event.status,
            "ORCAMON_MESSAGE": event.text,
        }
        try:
            self._children.append(subprocess.Popen(
                self.command, shell=True, env=env, stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True,
            ))
        except OSError:
            pass
