from trigger_mrna_generator import trigger_candidate

best_seq, best_score = trigger_candidate(36, "GGAAGGAGGUAACAAUG")
print("Best sequence:", best_seq, "Score:", best_score)