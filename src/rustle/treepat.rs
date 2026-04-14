//! Trie-like structure for transfrag lookup (C++ reference CTreePat).

use crate::graph::Graph;

#[derive(Debug)]
pub struct TreePat {
    _nodeno: usize,
    childno: usize,
    tr_idx: Option<usize>,
    next: Vec<Option<Box<TreePat>>>,
}

impl TreePat {
    pub fn new(nodeno: usize, childno: usize) -> Self {
        let mut next: Vec<Option<Box<TreePat>>> = Vec::with_capacity(childno);
        next.resize_with(childno, || None);
        Self {
            _nodeno: nodeno,
            childno,
            tr_idx: None,
            next,
        }
    }

    fn ensure_child(&mut self, idx: usize, nodeno: usize, childno: usize) -> Option<&mut TreePat> {
        if idx >= self.childno {
            return None;
        }
        if self.next[idx].is_none() {
            self.next[idx] = Some(Box::new(TreePat::new(nodeno, childno)));
        }
        self.next[idx].as_deref_mut()
    }

    fn get_child(&self, idx: usize) -> Option<&TreePat> {
        if idx >= self.childno {
            return None;
        }
        self.next[idx].as_deref()
    }
}

#[derive(Debug)]
pub struct TreePatIndex {
    root: TreePat,
    gno: usize,
}

impl TreePatIndex {
    pub fn new(gno: usize) -> Self {
        let root = TreePat::new(0, gno.saturating_sub(1));
        Self { root, gno }
    }

    fn child_index(prev: Option<usize>, next: usize, has_edge: bool, gno: usize) -> Option<usize> {
        if next == 0 {
            return None;
        }
        if let Some(p) = prev {
            if next <= p {
                return None;
            }
            let base = next - p - 1;
            if has_edge {
                Some((gno - 1 - p) + base)
            } else {
                Some(base)
            }
        } else {
            Some(next - 1)
        }
    }

    fn childno_for(next: usize, gno: usize) -> usize {
        2 * (gno.saturating_sub(next + 1))
    }

    pub fn find(
        &self,
        graph: &Graph,
        node_ids: &[usize],
        pattern: &crate::bitvec::GBitVec,
    ) -> Option<usize> {
        if node_ids.is_empty() {
            return None;
        }
        let mut tree = &self.root;
        let mut prev: Option<usize> = None;
        for &n in node_ids {
            let has_edge = prev
                .and_then(|p| graph.edge_bit_index(p, n))
                .map(|eid| pattern.get_bit(eid))
                .unwrap_or(false);
            let idx = Self::child_index(prev, n, has_edge, self.gno)?;
            tree = tree.get_child(idx)?;
            prev = Some(n);
        }
        tree.tr_idx
    }

    pub fn insert(
        &mut self,
        graph: &Graph,
        node_ids: &[usize],
        pattern: &crate::bitvec::GBitVec,
        tr_idx: usize,
    ) -> bool {
        if node_ids.is_empty() {
            return false;
        }
        let mut tree = &mut self.root;
        let mut prev: Option<usize> = None;
        for &n in node_ids {
            let has_edge = prev
                .and_then(|p| graph.edge_bit_index(p, n))
                .map(|eid| pattern.get_bit(eid))
                .unwrap_or(false);
            let idx = match Self::child_index(prev, n, has_edge, self.gno) {
                Some(v) => v,
                None => return false,
            };
            let childno = Self::childno_for(n, self.gno);
            let next = match tree.ensure_child(idx, n, childno) {
                Some(t) => t,
                None => return false,
            };
            tree = next;
            prev = Some(n);
        }
        tree.tr_idx = Some(tr_idx);
        true
    }
}
