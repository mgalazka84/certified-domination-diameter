# Certified domination in bipartite graphs

This repository reproduces the finite checks in *When does the certified domination number equal the size of the smaller part?* by Marek Gałązka and Hanna Wdowicka. A dominating set is certified if none of its vertices has exactly one neighbor outside the set. The paper asks when the minimum size of a certified dominating set equals the size of the smaller part of a bipartite graph.

## Run the checks

From the project root, use Python 3.12 and run

```sh
python -m pip install -r requirements.txt
python verification/verify.py
```

The output can be compared with [verification/results.txt](verification/results.txt). It was obtained with Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. The reduction check uses the HiGHS mixed-integer solver included in SciPy. To run a single check, pass one of the following arguments.

| Argument | Computation | Recorded result |
|---|---|---|
| `diameter` | Compute distances by breadth-first search and check certified dominating sets by subset enumeration for the fixed part sizes in the output. | Of 122,271 labeled graphs of diameter at most three, 174 satisfy the equality. All 106,866 graphs with part sizes (4,4) or (4,5) admit a smaller set. |
| `algorithm` | Compare the small-part dynamic program with exhaustive minimization for part sizes (1,3), (2,3), (3,3) and (3,4). | The minima agree on all 4,680 labeled graphs, including disconnected examples. |
| `reduction` | Enumerate isolate-free graphs on two through five vertices up to isomorphism. For each target size, compute the input vertex cover number by exhaustive search and the output certified domination number by integer programming. | All 184 instances from 33 input graphs satisfy the claimed relation and distance bounds. The one-edge example has 13 vertices and exactly two minimum sets, each of size six. |

The first two checks use only the Python standard library. They group graphs that differ only by a permutation of labels in the larger part and count them with their multiplicities.

## Verification code

OpenAI's ChatGPT assisted in writing the Python code used for the computational checks.
