/**
 * Author: Oleg
 * Date: 2024-11-13
 * License: CC0
 * Source: Oleg
 * Description: Gaussian elimination
 */

// Initialize augmented matrix
// Each row has variables +1 entries
vector<vector<int>> aug_matrix(equations, vector<int>(variables +1, 0));

// Populate augmented matrix
for(int i=0; i < k; ++i){
    for(int j=0; j < n; ++j){
        // Equation for C_blocks[i][j] = sum M[j][m] * P_blocks[i][m] mod 37
        for(int m=0; m < n; ++m){
            int var_index = j * n + m; // M[j][m]
            aug_matrix[i*n + j][var_index] = P_blocks[i][m];
        }
        // Set the augmented value
        aug_matrix[i*n + j][variables] = C_blocks[i][j];
    }
}

int rank =0;
vector<int> where(variables, -1);

for(int col=0; col < variables && rank < equations; ++col){
    // Find pivot
    int sel = -1;
    for(int i=rank; i < equations; ++i){
        if(aug_matrix[i][col] !=0){
            sel = i;
            break;
        }
    }
    if(sel == -1){
        continue; // No pivot in this column
    }
    swap(aug_matrix[sel], aug_matrix[rank]);

    int inv = mod_inverse(aug_matrix[rank][col], MOD);
    if(inv ==0){
        cout << "No solution\n";
        return 0;
    }
    for(int j=col; j <= variables; ++j){
        aug_matrix[rank][j] = (1LL * aug_matrix[rank][j] * inv) % MOD;
    }

    for(int i=0; i < equations; ++i){
        if(i != rank && aug_matrix[i][col] !=0){
            int factor = aug_matrix[i][col];
            for(int j=col; j <= variables; ++j){
                aug_matrix[i][j] = (aug_matrix[i][j] - 1LL * factor * aug_matrix[rank][j]) % MOD;
                if(aug_matrix[i][j] <0) aug_matrix[i][j] += MOD;
            }
        }
    }

    where[col] = rank;
    rank++;
}

// Check for inconsistency
for(int i=rank; i < equations; ++i){
    if(aug_matrix[i][variables] !=0){
        cout << "No solution\n";
        return 0;
    }
}

// Check if unique solution
bool unique = true;
for(int i=0; i < variables; ++i){
    if(where[i]==-1){
        unique = false;
        break;
    }
}

if(!unique){
    cout << "Too many solutions\n";
    return 0;
}

// Extract solution
vector<int> solution(variables, 0);
for(int i=0; i < variables; ++i){
    if(where[i]!=-1){
        solution[i] = aug_matrix[where[i]][variables];
    }
}

// Now, construct the matrix M from the solution
// M[j][m] where j = row, m = column
vector<vector<int>> M(n, vector<int>(n,0));
for(int j=0; j <n; ++j){
    for(int m=0; m <n; ++m){
        int var_index = j * n + m;
        M[j][m] = solution[var_index];
    }
}

// Output the matrix
for(int j=0; j <n; ++j){
    for(int m=0; m <n; ++m){
        cout << M[j][m];
        if(m != n-1){
            cout << ' ';
        }
    }
    cout << '\n';
}
