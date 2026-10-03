# Register FSANZ 2.9.4 Engine
try:
    fsanz_instance = FSANZ294Engine()
    registry.register(fsanz_instance)
except Exception as e:
    print(f"[Warning] FSANZ Engine registration error: {e}")

# Register Cutaneous Bioactivation Adapter
try:
    cutaneous_instance = CutaneousEngineAdapter()
    registry.register(cutaneous_instance)
except Exception as e:
    print(f"[Warning] Cutaneous Engine registration error: {e}")
    
