using BinaryFields, BinaryReedSolomon
using ZODATensor

n = 2^12
inv_rate = 4

# Data to be committed to
X_tilde = rand(BinaryElem16, n, n)

# Get the prover and verifier commitments
@info "--- Encoding+commitment time and commitment size"
@time prover_comm, verifier_comm = commit(X_tilde; inv_rate)

@show sizeof(verifier_comm)

# Get randomness (usually from Fiat-Shamir), also XXX don't be dumb
r = rand(BinaryElem16, n)
r_prime = r

@info "--- Encoding proof time and size"
@time enc_proof = encoding_proof(prover_comm, r, r_prime)
@show Base.format_bytes(sizeof(enc_proof))

# Sample some rows/cols
# S_row = sort(collect(Set(rand(1:n*inv_rate, 148))))
# S_col = sort(collect(Set(rand(1:n*inv_rate, 148))))
S_row = sort(collect(Set(rand(1:n*inv_rate, 24))))
S_col = sort(collect(Set(rand(1:n*inv_rate, 24))))

@info "--- Row/col openings size"
rowcol_proof = rowcol_openings(prover_comm, S_row, S_col)
@show Base.format_bytes(sizeof(rowcol_proof))

rs = reed_solomon(BinaryElem16, n, inv_rate*n)
ZODATensor.verify(verifier_comm, enc_proof, S_row, S_col, rowcol_proof, r, r_prime; rs)