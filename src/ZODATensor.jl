module ZODATensor

# Pseudocode: 
# prover commitment: X_tilde -> encode rows/cols -> X -> commit to rows/cols of X -> 2 merkle roots
# Noninteractive proof part? prover receive some randomness (fiat-shamir/could be anything) r and r_prime -> y_r = X_tilde*r, w_r_prime = X_tilde'*r_prime
# Interactive sampling? Verifier samples rows/cols of X, X_S are the S rows and Y_S' are the (supposed) columns of X
# Verifier check: X = G*X_tilde*G', _surely_ it is the case that X[:, 1:size(X_tilde)[1]]*r = G*X_tilde*r = Gy_r, let's look at a few (X_S*r) =?= (G*y_r)_S
# We will then do the same for "Y" (which is supposed to be the columns of X) -> verifier knows that prover has some matrix Y_tilde
# in their head.
# y_r' * r_prime = w_r_prime' * r
# If we do all of these check then => verifier is convinced that X_tilde exists, is consistent, and rows S/cols S' received are good

using BinaryFields, BinaryReedSolomon, BatchedMerkleTree

include("./zodaprover.jl")

end # module ZODATensor
