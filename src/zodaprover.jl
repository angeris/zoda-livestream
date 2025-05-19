export commit, encoding_proof, rowcol_openings

struct ZODAProverCommitment{T}
    X_tilde::Matrix{T}
    X::Matrix{T}
    merkle_row::CompleteMerkleTree
    merkle_col::CompleteMerkleTree
end

Base.show(io::IO, zoda::ZODAProverCommitment) = print(io, "ZODAProverCommitment{$(eltype(zoda.X_tilde))}")

struct ZODAVerifierCommitment
    merkle_row_root::MerkleRoot
    merkle_col_root::MerkleRoot
end

Base.sizeof(zoda::ZODAVerifierCommitment) = sizeof(zoda.merkle_row_root) + sizeof(zoda.merkle_col_root)
Base.show(io::IO, ::ZODAVerifierCommitment) = print(io, "ZODAVerifierCommitment")

# This thing is over BinaryElem32
function commit(X_tilde; inv_rate=4, rs_init=nothing)
    n, np = size(X_tilde)
    @assert n == np
    
    rs = isnothing(rs_init) ? reed_solomon(eltype(X_tilde), n, inv_rate*n) : rs_init

    X_cols = hcat(encode.(rs, eachcol(X_tilde))...)
    X = hcat(encode.(rs, eachrow(X_cols))...)'

    merkle_row = build_merkle_tree(eachrow(X))
    merkle_col = build_merkle_tree(eachcol(X))

    zoda_prover_commitment = ZODAProverCommitment{eltype(X_tilde)}(X_tilde, X, merkle_row, merkle_col)
    zoda_verifier_commitment = ZODAVerifierCommitment(get_root(merkle_row), get_root(merkle_col))

    return zoda_prover_commitment, zoda_verifier_commitment
end

struct ZODAEncodingProof{T}
    y_r::Vector{T}
    w_r_prime::Vector{T}
end

Base.sizeof(zoda::ZODAEncodingProof) = sizeof(zoda.y_r) + sizeof(zoda.w_r_prime)
Base.show(io::IO, zoda::ZODAEncodingProof{T}) where T = print(io, "ZODAEncodingProof{$(T)}")

# The randomness to be drawn from BinaryElem128
function encoding_proof(prover_commit::ZODAProverCommitment{T}, r, r_prime) where T
    ZODAEncodingProof(prover_commit.X_tilde * r, prover_commit.X_tilde' * r_prime)
end

struct ZODARowColOpenings{T}
    row_openings::Matrix{T}
    col_openings::Matrix{T}
    row_opening_proof::BatchedMerkleProof
    col_opening_proof::BatchedMerkleProof
end

Base.show(io::IO, zoda::ZODARowColOpenings) = print(io, "ZODARowColOpenings{$(eltype(zoda.row_openings))}")
Base.sizeof(zoda::ZODARowColOpenings) = sizeof(zoda.row_openings) + sizeof(zoda.col_openings) + sizeof(zoda.row_opening_proof) + sizeof(zoda.col_opening_proof)

function rowcol_openings(prover_commit::ZODAProverCommitment{T}, S_row, S_col) where T
    row_opening_proof = prove(prover_commit.merkle_row, S_row)
    col_opening_proof = prove(prover_commit.merkle_col, S_col)

    n = size(prover_commit.X_tilde, 1)

    row_openings = prover_commit.X[S_row, 1:n]
    col_openings = collect(prover_commit.X[1:n, S_col]')

    return ZODARowColOpenings{T}(row_openings, col_openings, row_opening_proof, col_opening_proof)
end

# Assume S and S_prime are in sorted order
function verify(verifier_commit::ZODAVerifierCommitment, verifier_enc_proof::ZODAEncodingProof{T}, S_row, S_col, rc_openings::ZODARowColOpenings{T}, r, r_prime; rs) where T
    y_r = verifier_enc_proof.y_r
    w_r_prime = verifier_enc_proof.w_r_prime

    # S_opening_paths, S_prime_opening_paths = merkle_opening_proofs
    X_S, Y_S_prime = rc_openings.row_openings, rc_openings.col_openings

    if X_S*r != encode(rs, y_r)[S_row]
        @error "X_S*r != encode(rs, y_r)[S]"
    end
    if Y_S_prime*r_prime != encode(rs, w_r_prime)[S_col]
        @error "Y_S_prime*r_prime != encode(rs, w_r_prime)[S]"
    end

    tree_depth = log_block_length(rs)

    # Check the rows
    if !BatchedMerkleTree.verify(
        verifier_commit.merkle_row_root,
        rc_openings.row_opening_proof;
        depth=tree_depth,
        leaves=encode.(rs, eachrow(X_S)),
        leaf_indices=S_row,
    )
        @error "Row opening proof failed"
    end

    # Check the rows
    if !BatchedMerkleTree.verify(
        verifier_commit.merkle_col_root,
        rc_openings.col_opening_proof;
        depth=tree_depth,
        leaves=encode.(rs, eachrow(Y_S_prime)),
        leaf_indices=S_col,
    )
        @error "Col opening proof failed"
    end

    if y_r' * r_prime != w_r_prime' * r
        @error "y_r'*r_prime != w_r_prime'*r"
    end
end