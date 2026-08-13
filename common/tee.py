#!/usr/bin/env python3
"""標準出力とログファイルへ同時に書き出す。"""


class Tee:
    """標準出力への出力を results.txt にも同時に書き出すためのクラス"""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()
