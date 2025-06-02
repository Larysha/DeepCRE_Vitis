# DeepCRE-Vitis

This repository contains an adapted version of the deepCRE re-implementation, customised specifically for *Vitis vinifera* (grapevine) to investigate gene regulatory mechanisms underlying water stress adaptation.

The original deepCRE project ([NAMlab/deepCRE_reimplemented](https://github.com/NAMlab/deepCRE_reimplemented)) is a deep learning framework designed to predict cis-regulatory elements and understand gene regulation in plants using sequence-based models.

---

This fork focuses on tailoring the deepCRE architecture and workflow to the *Vitis vinifera* genome and transcriptomic data. Key adaptations include:

- Integration of the *Vitis vinifera* T2T genome assembly and its specific gene annotations.
- Retraining the deepCRE model using grapevine-specific expression data to improve predictive performance in this species.
- Other possible adaptations down the line...

---

## Project Motivation

Water stress is a major environmental challenge affecting grapevine growth and yield. Understanding the cis-regulatory landscape that governs drought response genes is critical to developing resilient grapevine cultivars.