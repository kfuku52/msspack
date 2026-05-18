# Third-Party Notices

The `msspack` repository itself is distributed under the MIT License in [`LICENSE`](LICENSE). This file records attribution and preserved notices for third-party code or logic that was adapted into the current implementation.

## GFF2MSS

`msspack` includes logic in `src/msspack/mss_converter/` that was originally adapted from `GFF2MSS` and then refactored into native `msspack` modules. `msspack` is not distributed as a wrapper around the original `GFF2MSS` package, but attribution for the adapted converter logic is retained here.

Affected area:
- `src/msspack/mss_converter/`

Original project:
- [GFF2MSS](https://github.com/maedat/GFF2MSS)
- Copyright 2019-2020 Taro Maeda

License:

```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
