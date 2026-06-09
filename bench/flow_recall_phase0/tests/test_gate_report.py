import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gate_report as gr

def test_gate_proceed():
    v = gr.decide(recall_ceiling=200, best_auc=0.82)
    assert v["verdict"] == "PROCEED"

def test_gate_stop_no_discriminator():
    v = gr.decide(recall_ceiling=200, best_auc=0.55)
    assert v["verdict"] == "STOP"
    assert "discriminator" in v["reason"]

def test_gate_stop_small_ceiling():
    v = gr.decide(recall_ceiling=40, best_auc=0.9)
    assert v["verdict"] == "STOP"
    assert "ceiling" in v["reason"]
