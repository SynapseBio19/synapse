"""Synapse automated test suite (standard-library ``unittest``).

Run from the project root:

    python -m unittest discover -s tests -t .

Tests never contact the network: ``tests._support`` blocks sockets.
"""
