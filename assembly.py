def k_mers(k, *substrings):
    kmers = []
    for substring in substrings:
        for i in range(len(substring)):
            curr_substring = substring[i:i+k]
            if len(curr_substring) == k:
                kmers.append(curr_substring)
            else:
                break
    return kmers

def prefixes_suffixes(kmers):
    prefs = [k[:-1] for k in kmers]
    suffs = [k[1:] for k in kmers]
    return prefs + suffs
        
