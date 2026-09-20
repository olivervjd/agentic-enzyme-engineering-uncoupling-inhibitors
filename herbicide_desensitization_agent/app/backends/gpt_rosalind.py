"""Backward-compatible imports; live adapters use the configured model, now GPT-5.6 Luna."""
from .gpt_json import GPTJSONBackend

GPTRosalindJSONBackend = GPTJSONBackend
