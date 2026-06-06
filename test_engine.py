"""
Test script for IntelliPrime engine.
Run this to verify everything works.
"""

from engine import main

print("Starting IntelliPrime test...\n")

# Test with Human Hemoglobin Beta (HBB)
result = main("NM_000518.5")

# Verify the output
print("\n--- Quick Validation ---")
assert result["cds_length"] > 0, "ERROR: CDS not extracted"
assert len(result["optimized_cds"]) == result["cds_length"], "ERROR: Codon optimization changed length"
assert "CATATG" in result["complete_template"], "ERROR: NdeI site missing"
assert "CTCGAG" in result["complete_template"], "ERROR: XhoI site missing"

print("✅ All assertions passed!")
print("✅ Engine is working correctly.")