"""Tests for read-to-gene assignment.

The policy under test is the one the predecessor got wrong twice: once by
dropping ties between genes that are not separable, and once by discarding the
record of which pairs those were. Both failures are silent — the output looks
like a gene with poor coverage rather than like a bug — so the tests assert the
behaviour explicitly rather than trusting it to surface downstream.
"""

from __future__ import annotations

import pytest

from lilrwgs.assign import (
    DEFAULT_SHARED_GROUPS,
    arbitrate,
    group_index,
    read_shared_names,
    write_shared_table,
)


class TestGroupIndex:
    def test_finds_a_member(self):
        assert group_index("LILRA6", DEFAULT_SHARED_GROUPS) == 0
        assert group_index("LILRB1", DEFAULT_SHARED_GROUPS) == 1

    def test_separable_gene_has_no_group(self):
        assert group_index("LILRA4", DEFAULT_SHARED_GROUPS) is None


class TestArbitrate:
    def test_single_hit_is_assigned(self):
        assigned, stats = arbitrate({"LILRA4": {"r1": -10}, "LILRA5": {}})
        assert assigned["LILRA4"]["r1"] == 1
        assert stats.n_unique == 1

    def test_clear_winner_takes_the_pair(self):
        assigned, stats = arbitrate({
            "LILRA4": {"r1": -2},
            "LILRA5": {"r1": -30},
        })
        assert assigned["LILRA4"]["r1"] == 1
        assert "r1" not in assigned["LILRA5"]
        assert stats.n_arbitrated == 1

    def test_tie_across_separable_genes_is_dropped(self):
        """Where the references are separable, a tie means the read genuinely
        cannot be placed, and keeping it in both would invent coverage."""
        assigned, stats = arbitrate({
            "LILRA4": {"r1": -2},
            "LILRA5": {"r1": -2},
        })
        assert not assigned["LILRA4"] and not assigned["LILRA5"]
        assert stats.n_dropped == 1

    def test_tie_within_a_shared_group_is_kept_in_both(self):
        """The central correction. LILRA6/LILRB3 share a ~4.6 kb block with no
        gene-diagnostic 31-mers; a competitive discard deletes it from both
        genes at once, which cost 26% of the LILRB3 CDS in the capture cohort.
        """
        assigned, stats = arbitrate({
            "LILRA6": {"r1": -2},
            "LILRB3": {"r1": -2},
        })
        assert assigned["LILRA6"]["r1"] == 2
        assert assigned["LILRB3"]["r1"] == 2
        assert stats.n_shared == 1
        assert stats.n_dropped == 0

    def test_tie_spanning_two_groups_is_still_dropped(self):
        """Being in *a* shared group is not enough — the tie has to be confined
        to one, or it is genuinely ambiguous."""
        assigned, stats = arbitrate({
            "LILRA6": {"r1": -2},
            "LILRB1": {"r1": -2},
        })
        assert stats.n_dropped == 1
        assert not assigned["LILRA6"]

    def test_tie_between_a_grouped_and_ungrouped_gene_is_dropped(self):
        assigned, stats = arbitrate({
            "LILRA6": {"r1": -2},
            "LILRA4": {"r1": -2},
        })
        assert stats.n_dropped == 1

    def test_buffer_controls_what_counts_as_tied(self):
        scores = {"LILRA6": {"r1": -2}, "LILRB3": {"r1": -5}}
        wide, wide_stats = arbitrate(scores, buffer=5)
        assert wide_stats.n_shared == 1
        narrow, narrow_stats = arbitrate(scores, buffer=1)
        assert narrow_stats.n_arbitrated == 1
        assert narrow["LILRA6"]["r1"] == 1

    def test_three_way_tie_inside_one_group_keeps_all(self):
        groups = [("A", "B", "C")]
        assigned, stats = arbitrate(
            {"A": {"r1": -2}, "B": {"r1": -2}, "C": {"r1": -2}},
            shared_groups=groups,
        )
        assert all(assigned[g]["r1"] == 3 for g in "ABC")
        assert stats.n_shared == 1

    def test_shared_groups_can_be_disabled(self):
        """Passing an empty list must mean 'no shared groups', not 'use the
        defaults' — otherwise the policy cannot be turned off for a comparison.
        """
        _, stats = arbitrate({"LILRA6": {"r1": -2}, "LILRB3": {"r1": -2}},
                             shared_groups=[])
        assert stats.n_dropped == 1
        assert stats.n_shared == 0


class TestStats:
    def test_counts_every_pair_exactly_once(self):
        assigned, stats = arbitrate({
            "LILRA4": {"solo": -3, "win": -2, "sep": -2},
            "LILRA5": {"win": -30, "sep": -2},
            "LILRA6": {"share": -2},
            "LILRB3": {"share": -2},
        })
        assert stats.n_total == 4
        assert (stats.n_unique, stats.n_arbitrated,
                stats.n_shared, stats.n_dropped) == (1, 1, 1, 1)

    def test_drop_rate_is_visible(self):
        """The predecessor's driver piped the arbitration log through
        `grep -vF "Ambiguous read pair"`, hiding that 70.9% of LILRA6 pairs were
        being discarded. The number has to be in the output, not the log."""
        _, stats = arbitrate({
            "LILRA4": {f"r{i}": -2 for i in range(9)},
            "LILRA5": {f"r{i}": -2 for i in range(9)},
        })
        assert stats.as_row()["dropped_fraction"] == 1.0

    def test_per_gene_counts_include_shared(self):
        _, stats = arbitrate({"LILRA6": {"r1": -2}, "LILRB3": {"r1": -2}})
        assert stats.per_gene["LILRA6"] == 1
        assert stats.per_gene_shared["LILRA6"] == 1


class TestSharedTable:
    def test_round_trips_through_the_fastq_hop(self, tmp_path):
        """Reads are recruited against a pangenome panel and called against a
        single reference, so there is a FASTQ in between and BAM tags do not
        survive it. This table is how the shared-block fact crosses that gap —
        without it the depth model cannot tell a legitimate two-gene block from
        an unseparated pile-up.
        """
        assigned, _ = arbitrate({
            "LILRA6": {"shared": -2, "own": -1},
            "LILRB3": {"shared": -2},
        })
        path = tmp_path / "shared.tsv"
        n = write_shared_table(assigned, path)
        assert n == 2

        a6 = read_shared_names(path, "LILRA6")
        assert a6 == {"shared": 2}
        assert "own" not in a6, "uniquely assigned pairs must not be marked shared"

    def test_missing_table_is_empty_not_an_error(self, tmp_path):
        assert read_shared_names(tmp_path / "absent.tsv", "LILRA6") == {}
