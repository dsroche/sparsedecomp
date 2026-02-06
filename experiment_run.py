import random
from pathlib import Path
from typing import Self
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import cattrs
from sklearn.metrics import roc_curve
from uniplot import plot
import plotext
from sage.all import (
    GF,
    PolynomialRing,
    next_prime,
)
from progbar import ProgressBar
from dense_decomp import (
    collisions,
    rand_poly,
    rand_decomp,
    rand_indecomp,
)

# Register how to turn an array into a list (for saving)
cattrs.register_unstructure_hook(
    np.ndarray, lambda v: v.tolist()
)

# Register how to turn a list back into an array (for loading)
cattrs.register_structure_hook(
    np.ndarray, lambda v, _: np.array(v)
)

@dataclass
class Range:
    low: int
    high: int

    def sample(self) -> int:
        return random.randrange(self.low, self.high+1)

@dataclass
class Poly:
    coeffs: list[int]

    def to_sage(self, PR):
        return PR(coeffs)

    @classmethod
    def from_sage(cls, f) -> Self:
        return cls(coeffs=[int(c) for c in f.list()])

@dataclass
class DataCase:
    poly: Poly
    collisions: list[int]

@dataclass
class DataSet:
    decomp: list[DataCase]
    indecomp: list[DataCase]

@dataclass
class Experiment:
    p: int
    rrange: Range
    drange: Range
    trange: Range|None
    nsamples: int
    ntrain: int
    ntest: int

    train_data: DataSet|None = None
    test_data: DataSet|None = None

    thresholds: list[float]|None = None
    accuracies: list[float]|None = None

    def same_params(self, other: 'Experiment') -> bool:
        return (
            (self.p, self.rrange, self.drange, self.trange, self.nsamples, self.ntrain, self.ntest)
            ==
            (other.p, other.rrange, other.drange, other.trange, other.nsamples, other.ntrain, other.ntest)
            )

    @cached_property
    def F(self):
        return GF(self.p)

    @cached_property
    def PR(self):
        return PolynomialRing(self.F, 'x')

    def build_one(self) -> tuple[DataCase, DataCase]:
        while True:
            r = self.rrange.sample()
            s = round(self.drange.sample() / r)
            if r < 2 or s < 2:
                continue
            break
        gt = r+1 if self.trange is None else self.trange.sample()
        g = rand_poly(self.PR, r, gt)
        ht = s+1 if self.trange is None else self.trange.sample()
        h = rand_poly(self.PR, s, ht)
        decom = g(h)
        decol = collisions(decom, self.F, self.nsamples)

        indecom = rand_indecomp(self.PR, decom.degree(), decom.number_of_terms())
        indecol = collisions(indecom, self.F, self.nsamples)

        return (DataCase(Poly.from_sage(decom), decol),
                DataCase(Poly.from_sage(indecom), indecol))

    def build_dataset(self, n: int, ppe: ProcessPoolExecutor, pb: ProgressBar) -> DataSet:
        dcases = []
        icases = []
        futures = [ppe.submit(self.build_one) for _ in range(n)]
        for fut in as_completed(futures):
            dcase, icase = fut.result()
            dcases.append(dcase)
            icases.append(icase)
            pb += 1
        return DataSet(dcases, icases)

    def build_data(self) -> None:
        tot = self.ntrain + self.ntest
        print(f"Building {tot} test cases...")
        with ProgressBar(tot) as pb:
            with ProcessPoolExecutor() as ppe:
                self.train_data = self.build_dataset(self.ntrain, ppe, pb)
                self.test_data = self.build_dataset(self.ntest, ppe, pb)

    def train(self) -> None:
        assert self.train_data is not None
        print("Training...")
        with ProgressBar(self.nsamples+1) as pb:
            self.thresholds = []
            decomps = np.array([dc.collisions for dc in self.train_data.decomp]).T
            indecomps = np.array([dc.collisions for dc in self.train_data.indecomp]).T
            yt = np.concat([np.zeros(self.ntrain), np.ones(self.ntrain)])
            for dec, indec in zip(decomps, indecomps):
                fpr,tpr,thresh = roc_curve(yt, np.concat([indec, dec]))
                ot = thresh[np.argmax(tpr-fpr)]
                try:
                    below = indec[indec < ot].max()
                    above = indec[indec >= ot].min()
                    cutoff = (below + above) / 2
                except ValueError:
                    cutoff = ot
                self.thresholds.append(float(cutoff))
                pb += 1

    def test(self) -> None:
        assert self.test_data is not None
        assert self.thresholds is not None
        print("Testing...")
        decomps = np.array([dc.collisions for dc in self.test_data.decomp]).T
        indecomps = np.array([dc.collisions for dc in self.test_data.indecomp]).T
        self.accuracies = []
        for cutoff, dec, indec in zip(self.thresholds, decomps, indecomps):
            true_pos = sum(1 for x in dec if x >= cutoff)
            true_neg = sum(1 for x in indec if x < cutoff)
            acc = (true_pos + true_neg) / (2 * self.ntest)
            self.accuracies.append(acc)

    def build_all(self, name: str) -> None:
        if self.accuracies:
            return
        print(f"Building {name}...")
        if not self.train_data or not self.test_data:
            self.build_data()
        if not self.thresholds:
            self.train()
        if not self.accuracies:
            self.test()
        print(f"done building {name}")

if __name__ == '__main__':
    expers = {
        'sparse_small': Experiment(
            p = int(next_prime(2**10)),
            rrange = Range(2, 16),
            drange = Range(500,1000),
            trange = Range(3, 5),
            nsamples = 300,
            ntrain = 1000,
            ntest =  1000,
        ),
        'dense_small': Experiment(
            p = int(next_prime(2**10)),
            rrange = Range(2, 200),
            drange = Range(500,1000),
            trange = None,
            nsamples = 300,
            ntrain = 1000,
            ntest =  1000,
        ),
        'sparse_big': Experiment(
            p = int(next_prime(2**16)),
            rrange = Range(2, 16),
            drange = Range(500,1000),
            trange = Range(3, 5),
            nsamples = 2400,
            ntrain = 500,
            ntest =  500,
        ),
        'dense_big': Experiment(
            p = int(next_prime(2**16)),
            rrange = Range(2, 200),
            drange = Range(500,1000),
            trange = None,
            nsamples = 2400,
            ntrain = 500,
            ntest =  500,
        ),
    }
    data = Path('data')
    data.mkdir(exist_ok=True)
    for name, exp in expers.items():
        f = data / f"{name}.json"
        if f.exists():
            try:
                loaded = cattrs.structure(json.load(f.open()), Experiment)
                if exp.same_params(loaded):
                    expers[name] = loaded
                    print(f"loaded {name} from {f}")
                    continue
            except:
                pass
        exp.build_all(name)
        with f.open('w') as fout:
            json.dump(cattrs.unstructure(exp), fout)
        print(f"saved {name} to {f}")
    plotext.plotsize(plotext.terminal_width(), plotext.terminal_height()*.5)

    plotext.plot(expers['sparse_small'].accuracies, label="sparse")
    plotext.plot(expers['dense_small'].accuracies, label="dense")
    plotext.title("Small")
    plotext.show()

    plotext.clear_data()
    plotext.plot(expers['sparse_big'].accuracies, label="sparse")
    plotext.plot(expers['dense_big'].accuracies, label="dense")
    plotext.title("Big")
    plotext.show()
