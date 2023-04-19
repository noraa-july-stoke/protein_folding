import heapq
import json

with open("amino_acids.json") as f:
    amino_acids_data = json.load(f)

amino_acids = amino_acids_data["amino_acids"]

def calculate_score(bonded_amino_acids):
    score = 0
    for amino_acid in bonded_amino_acids:
        score += amino_acid["properties"]["polarity"] == "polar"
    return score

def assemble_protein(start_amino_acid, bonded_amino_acids, ratios):
    best_score = -float('inf')
    results = []
    heap = [(0, start_amino_acid, bonded_amino_acids, ratios)]

    while heap:
        score, amino_acid, bonded_amino_acids, ratios = heapq.heappop(heap)

        if len(bonded_amino_acids) == len(amino_acids):
            new_score = calculate_score(bonded_amino_acids)
            if new_score > best_score:
                results.append((bonded_amino_acids, new_score))
                best_score = new_score
        else:
            for bond_type in amino_acid["allowable_bond_types"]:
                for next_amino_acid in amino_acids:
                    if next_amino_acid != amino_acid and bond_type in next_amino_acid["allowable_bond_types"]:
                        if next_amino_acid not in bonded_amino_acids:
                            if check_ratio(next_amino_acid, ratios):
                                new_bonded_amino_acids = bonded_amino_acids + [next_amino_acid]
                                new_ratios = update_ratios(next_amino_acid, ratios)
                                new_score = calculate_score(new_bonded_amino_acids)
                                lower_bound = new_score - len(amino_acids) + len(new_bonded_amino_acids)
                                print("===================LOWER_BOUND===================")
                                print(-lower_bound)
                                print("===================LOWER_BOUND===================")

                                print("============NEXT AMINO ACID============")
                                print(next_amino_acid)
                                print("============NEXT AMINO ACID============")

                                print("============NEW BONDED AMINO ACIDS==========")
                                print(new_bonded_amino_acids)
                                print("============NEW BONDED AMINO ACIDS==========")

                                print("===========NEW RATIOS============")
                                print(new_ratios)
                                print("===========NEW RATIOS============")

                                heapq.heappush(heap, (-lower_bound, next_amino_acid, new_bonded_amino_acids, new_ratios))
                                print("HEAPQ PUSHED")

    return results


def calculate_lower_bound(ratios):
    max_score = 0
    for amino_acid in amino_acids:
        if check_ratio(amino_acid, ratios):
            max_score += calculate_score([amino_acid])
    return max_score


def calculate_score(bonded_amino_acids):
    total_score = 0
    for amino_acid in bonded_amino_acids:
        amino_acid_score = 0
        if amino_acid["properties"]["polarity"] == "polar":
            amino_acid_score += 1
        elif amino_acid["properties"]["polarity"] == "nonpolar":
            amino_acid_score -= 1
        if amino_acid["properties"]["charge"] == "positive":
            amino_acid_score += 1
        elif amino_acid["properties"]["charge"] == "negative":
            amino_acid_score -= 1
        total_score += amino_acid_score
    print(total_score)
    return total_score


def check_ratio(amino_acid, ratios):
    for key, value in ratios.items():
        if amino_acid["abbreviation"] == key:
            return value > 0
    return False


def update_ratios(amino_acid, ratios):
    updated_ratios = ratios.copy()
    for key, value in updated_ratios.items():
        if amino_acid["abbreviation"] == key:
            updated_ratios[key] = value - 1
    return updated_ratios


# define the ratios of amino acids
ratios = {
    "Ala": 5,
    "Gly": 3,
    "Ser": 2,
    "Val": 2,
    "Asn": 1,
    "Asp": 1,
    "Gln": 1,
    "Glu": 1,
    "His": 1,
    "Ile": 1,
    "Leu": 1,
    "Lys": 1,
    "Met": 1,
    "Phe": 1,
    "Pro": 1,
    "Thr": 1,
    "Trp": 1,
    "Tyr": 1,
    "Cys": 1
}

start_amino_acid = amino_acids[0]
bonded_amino_acids = [start_amino_acid]
results = assemble_protein(start_amino_acid, bonded_amino_acids, ratios)
# print(results)
counter = 0

for result in results:
    # print("Bonded amino acids:", [amino_acid["name"] for amino_acid in result[0]])
    print("Score:", result[1])
    counter += 1
    if counter == 5:
        break
