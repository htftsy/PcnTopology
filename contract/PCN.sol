// SPDX-License-Identifier: GPL-3.0

/* For research purpose only */
pragma solidity =0.8.27;

contract PCN { 

    uint32 public N;
    uint32 public P;
    uint32 public Q;
    uint32[][] public adj;
    address payable[] public addr;

    uint32 public iter;
    uint32 public chPerIndex;
    uint32 public custimizedchPerIndex;
    uint32 public dutychPerIndex;
    uint256 public amountPerCh;
    uint32[] public remainedCustimizedChs;

    struct matrice {
        uint32 a;
        uint32 b;
        uint32 c;
        uint32 d;
    }

    struct chState {
        uint32 u;
        uint32 v;
        uint32 k;
        uint256 R;
        int256 frz;
        uint256 I;
        bytes sigu;
        bytes sigv;
    }

    matrice[] public PGLMembers;
    matrice[] public SMembers;
    mapping(uint32 => uint32) public matriceToIndex; 
    mapping(uint32 => uint32) public iterNumToIndex;
    mapping(uint32 => mapping(uint32 => uint32)) public adjacentIndex;
    mapping(uint32 => mapping(uint32 => uint256)) public capacity;
    mapping(uint32 => mapping(uint32 => uint256)) public heightForWithdl;
    mapping(uint32 => mapping(uint32 => chState)) public stateForWithdl;

    function inverseinFq(uint32 v) 
        public view returns (uint32) 
    {
        // Since Q is small, it is faster to enumerate directly.
        for(uint32 i = 0; i < Q; i++)
            if(v * i % Q == 1)
                return i;
        return 0;
    }

    function matriceConstruct(uint32 a, uint32 b, uint32 c, uint32 d) 
        public view returns (matrice memory) 
    {
        matrice memory t;
        t.a = a % Q;
        t.b = b % Q;
        t.c = c % Q;
        t.d = d % Q;
        return t;
    }

    function nonZeroDeterminant(matrice memory x)
        public view returns (bool)
    {
        return (x.a * x.d + Q * Q - x.b * x.c) % Q != 0;
    }

    function normalize(matrice memory x)
        public view returns (matrice memory)
    {
        uint32 mulfac;
        if(x.a == 0) {
            if(x.b == 0) {
                if(x.c == 0) {
                    if(x.d == 0)
                        return x;
                    else
                        mulfac = inverseinFq(x.b);
                }
                else
                    mulfac = inverseinFq(x.c);
            }
            else
                mulfac = inverseinFq(x.b);
        }
        else
            mulfac = inverseinFq(x.a);
        x.a = x.a * mulfac % Q;
        x.b = x.b * mulfac % Q;
        x.c = x.c * mulfac % Q;
        x.d = x.d * mulfac % Q;
        return x;
    }

    function matriceMul(matrice memory x, matrice memory y) 
        public view returns (matrice memory)
    {
        matrice memory t;
        t.a = (x.a * y.a + x.b * y.c) % Q;
        t.b = (x.a * y.b + x.b * y.d) % Q;
        t.c = (x.c * y.a + x.d * y.c) % Q;
        t.d = (x.c * y.b + x.d * y.d) % Q;
        return normalize(t);
    }

    function serializeMatriceToInt(matrice memory x)
        public view returns (uint32) 
    {
        uint32 res = 0;
        res = x.a;
        res = res * Q + x.b;
        res = res * Q + x.c;
        res = res * Q + x.d;
        return res;
    }

    function genPGLMembers() public {
        matrice memory tmp;
        matrice memory tmp_normalized;
        uint32 ct = 0;
        for(uint32 i = 0; i < Q; i++)
            for(uint32 j = 0; j < Q; j++)
                for(uint32 k = 0; k < Q; k++)
                    for(uint32 l = 0; l < Q; l++) {
                        tmp.a = i; 
                        tmp.b = j;
                        tmp.c = k;
                        tmp.d = l;
                        if(! nonZeroDeterminant(tmp))
                            continue;
                        tmp_normalized = normalize(tmp);
                        if(tmp_normalized.a == i && tmp_normalized.b == j && tmp_normalized.c == k && tmp_normalized.d == l) {
                            PGLMembers.push(tmp);
                            ct ++;
                        }
                    }
        for(uint32 i = 0; i < N; i++)
            matriceToIndex[serializeMatriceToInt(PGLMembers[i])] = i;
    }

    function genSMembers() public {
        int32 x;
        int32 y;
        bool foundxy = false;
        int32 Pi32 = int32(P);
        int32 Qi32 = int32(Q);

        for(x = 0; x < Qi32; x++) {
            for(y = 0; y < Qi32; y++) {
                if((x * x + y * y + 1) % Qi32 == 0)
                    foundxy = true;
                if(foundxy) break;
            }
            if(foundxy) break;
        }
        require((x * x + y * y + 1) % Qi32 == 0);
        matrice memory t;
        uint32 ct = 0;
        for(int32 a0 = 0; a0 < Pi32; a0++)
            for(int32 a1 = 1 - Pi32; a1 < Pi32; a1++)
                for(int32 a2 = 1 - Pi32; a2 < Pi32; a2++)
                    for(int32 a3 = 1 - Pi32; a3 < Pi32; a3++) {
                        if(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3 != Pi32)
                            continue;
                        if(P % 4 == 1) {
                            if((a0 % 2) == 0 || a0 == 0)
                                continue;
                        }
                        else {
                            if(!((a0 % 2) == 0 && a0 > 0) && !(a0 == 0 && a1 > 0))
                                continue;
                        }
                        t.a = uint32(3 * Qi32 * Qi32 + a0 + a1 * x + a3 * y) % Q;
                        t.b = uint32(3 * Qi32 * Qi32 - a1 * y + a2 + a3 * x) % Q;
                        t.c = uint32(3 * Qi32 * Qi32 - a1 * y - a2 + a3 * x) % Q;
                        t.d = uint32(3 * Qi32 * Qi32 + a0 - a1 * x - a3 * y) % Q;
                        SMembers.push(normalize(t));
                        ct ++;
                    }
    }

    uint32 public front;
    uint32 public rear;
    uint32[] public qu;
    bool[] public vis;

    function setupLPS(uint32 i) public {

        for(uint32 j = 0; j <= uint32(P); j++) {
            matrice memory t;
            t = matriceMul(PGLMembers[i], SMembers[j]);
            uint32 k = matriceToIndex[serializeMatriceToInt(t)];
            adjacentIndex[i][j] = k;
        }
    }

    constructor (
        uint32 _P, 
        uint32 _Q, 
        uint32 custimized_Ch_per_index, // this is (1 - gamma) * m in the paper
        uint256 amount_per_Ch
    )  
    {
        uint32 Ch_per_index = _P + 1;
        require(custimized_Ch_per_index <= Ch_per_index);

        N = (_Q ** 3) - _Q;    // `**' means power for Solidity

        iter = 0;
        adj = new uint32[][](N);  
        addr = new address payable[](N);

        amountPerCh = amount_per_Ch;
        chPerIndex = Ch_per_index;
        custimizedchPerIndex = custimized_Ch_per_index;
        dutychPerIndex = chPerIndex - custimizedchPerIndex;

        remainedCustimizedChs = new uint32[](N);

        P = _P;
        Q = _Q;

        iterNumToIndex[0] = 0;
        front = 0;
        rear = 0;
        qu = new uint32[] (N + 1);
        qu[0] = 0;
        vis = new bool[] (N + 1);
        vis[0] = true;
        for(uint32 i = 1; i < N; i++)
            vis[i] = false;

        // To setup, please execute genPGLMembers() and genSMembers()
    }

    function corrdinateToIndex (matrice memory coor) 
        public 
        view 
        returns (uint32) 
    {
        return matriceToIndex[serializeMatriceToInt(coor)];
    }

    function indexToCoordinate (uint32 index) 
        public 
        view 
        returns (matrice memory coor)
    {
        return PGLMembers[index];
    }

    function registerIn ()
        public 
        payable 
        returns (uint32, uint32[] memory)
    {
        require(msg.value >= amountPerCh * uint256(chPerIndex));   // deposits for both duty Chs and typical Chs
        require(front < N);
        setupLPS(qu[front]);
        uint32 target;
        uint32 f = qu[front];
        for(uint32 j = 0; j <= P; j++) {
            target = adjacentIndex[f][j];
            if(! vis[target]) {
                qu[++ rear] = target;
                vis[target] = true;
                iterNumToIndex[rear] = target;
            }
        }
        front ++;

        uint32 index = iterNumToIndex[iter];
        iter ++;
        require(corrdinateToIndex(indexToCoordinate(index)) == index);

        addr[index] = payable(msg.sender);

        adj[index] = new uint32[](dutychPerIndex);

        for(uint32 k = 0; k < dutychPerIndex; k++)
                adj[index][k] = adjacentIndex[index][k];

        for(uint32 i = 0; i < dutychPerIndex; i++)
            capacity[index][uint32(adj[index][uint32(i)])] = amountPerCh;

        remainedCustimizedChs[index] = custimizedchPerIndex;

        return (index, adj[index]);
    }

    function addCustimizedCh (uint32 index, uint32 indexTo, bytes memory signature) 
        public 
    {
        require(addr[index] == msg.sender);
        require(remainedCustimizedChs[index] > 0);
        require(isValidSignature(addr[indexTo], uint256(index), signature));
        adj[index].push(indexTo);
        capacity[index][indexTo] = amountPerCh;
        remainedCustimizedChs[index] --;
    }

    function queryIndexAddress (uint32 index) 
        public 
        view 
        returns (address)
    {
        require(adj[index].length > 0); 
        return addr[index];
    }

    function queryAdjacentIndices (uint32 index)
        public
        view 
        returns (uint32[] memory) 
    {
        require(adj[index].length > 0); 
        return adj[index];
    }

    function showCapacity (uint32 indexSender, uint32 indexReceiver) 
        public 
        view 
        returns (uint256) 
    {
        return capacity[indexSender][indexReceiver];
    }

    function topUp (uint32 _indexSender, uint32 _indexReceiver, uint256 value) 
        public 
        payable 
    {
        uint32 indexSender = _indexSender;
        uint32 indexReceiver = _indexReceiver;
        require(addr[indexSender] == msg.sender);
        require(value <= msg.value);
        capacity[indexSender][indexReceiver] += value;
    }


    function symmetricChState(chState memory s) 
        internal 
        view
        returns (chState memory t)
    {
        t.u = s.v;
        t.v = s.u;
        t.k = s.k;
        t.R = 2 * amountPerCh - s.R;
        t.frz = - s.frz;
        t.I = s.I;
        t.sigu = s.sigv;
        t.sigv = t.sigu;
    }

    function hashChStateToSign(chState memory s)  
    // This function should evaluate to the same value for s and symmetricChState(s)
        internal 
        view 
        returns (uint256)
    {
        if(s.u > s.v)
            s = symmetricChState(s);
        return uint256(keccak256(abi.encode(s.u, s.v, s.k, s.R, s.frz, s.I))); 
    }

    function checkChStateValidity(chState memory s) 
        internal 
        view 
        returns (bool)
    {
        return capacity[s.u][s.v] > 0 && s.frz + int256(s.R) == int256(amountPerCh);
    }

    function verifyStateTransition(chState memory s, chState memory t) 
        internal 
        view 
        returns (bool)
    {
        if(s.u > s.v)
            s = symmetricChState(s);
        if(t.u > t.v)
            t = symmetricChState(t);
        if(s.u != t.u || s.v != t.v || t.k != s.k + 1 || s.frz > 0 && t.frz > 0)
            return false;
        if(t.sigu[0] > 0 && t.sigv[0] > 0) 
            return isValidSignature(addr[t.u], hashChStateToSign(t), t.sigu) 
                && isValidSignature(addr[t.v], hashChStateToSign(t), t.sigv);
        // Now it is unfreezing
        return s.R == t.R && s.frz != 0 && t.frz == 0 && uint256(keccak256(abi.encode(t.I))) == s.I
            && ( 
                (s.frz > 0 && isValidSignature(addr[t.v], hashChStateToSign(t), t.sigv)) 
             || (s.frz < 0 && isValidSignature(addr[t.u], hashChStateToSign(t), t.sigu))
            );
    }

    function withdlIssue(uint32 u, uint32 v, chState memory state_p, chState memory state) 
        public 
    {
        require(addr[u] == msg.sender);
        require(state.u == u && state.v == v);
        require(capacity[v][u] > 0);
        require(checkChStateValidity(state_p));
        require(checkChStateValidity(state));
        require(verifyStateTransition(state_p, state));
        require(heightForWithdl[u][v] == 0);

        heightForWithdl[u][v] = heightForWithdl[v][u] = block.number + 30;
        stateForWithdl[u][v] = stateForWithdl[v][u] = state;
    }

    function withdlComplete(uint32 u, uint32 v)
        public 
    {
        require(heightForWithdl[u][v] > 0);
        require(heightForWithdl[u][v] < block.number);

        bool found = false;
        uint32 len = uint32(adj[u].length);
        for(uint i = 0; i < len; i++) {
            if(uint32(adj[u][i]) == v) {
                found = true;
                adj[u][i] = adj[u][len - 1];
                adj[u].pop(); 
            }
        }
        require(found); 
        len = uint32(adj[v].length);
        for(uint i = 0; i < len; i++) {
            if(uint32(adj[v][i]) == u) {
                found = true;
                adj[v][i] = adj[v][len - 1];
                adj[v].pop(); 
            }
        }
        require(found); 
        uint256 amount = stateForWithdl[u][v].R;
        if(stateForWithdl[u][v].frz != 0)
            amount = uint256(int256(amount) + stateForWithdl[u][v].frz);
        if(stateForWithdl[u][v].u == u) {
            addr[u].transfer(amount);
            addr[v].transfer(2 * amountPerCh - amount);
        }
        else {
            addr[v].transfer(amount);
            addr[u].transfer(2 * amountPerCh - amount);
        }
        capacity[u][v] = capacity[v][u] = 0;
    }

    function disputeConflictRecord(uint32 u, uint32 v, chState memory witness) 
        public 
    {
        if(witness.v != v) 
            witness = symmetricChState(witness);
        require(addr[v] == msg.sender);
        require(witness.u == u && witness.v == v);
        require(witness.sigu[0] > 0 && witness.sigv[0] > 0);
        require(isValidSignature(addr[witness.u], hashChStateToSign(witness), witness.sigu));
        require(isValidSignature(addr[witness.v], hashChStateToSign(witness), witness.sigv));
        chState memory os = stateForWithdl[u][v];
        if(os.v != v)
            os = symmetricChState(os);
        require(os.u == u && os.v == v);
        if(witness.k == os.k) {
            if(witness.frz != 0 && os.frz == 0)
                return;
            stateForWithdl[u][v] = stateForWithdl[v][u] = witness;
        }
        else if (witness.k == os.k - 1) {
            if(witness.frz == 0 && os.sigv[0] == 0)
                stateForWithdl[u][v] = stateForWithdl[v][u] = witness;
        }
        else 
            require(false); // should not reach here
    }

    function disputeStateRenewal(uint32 u, uint32 v, chState memory ps, chState memory s) 
        public 
    {
        if(s.v != v) 
            s = symmetricChState(s);
        require(addr[v] == msg.sender);
        require(s.u == u && s.v == v);
        require(ps.sigu[0] > 0 && ps.sigv[0] > 0);
        require(isValidSignature(addr[ps.u], hashChStateToSign(ps), ps.sigu));
        require(isValidSignature(addr[ps.v], hashChStateToSign(ps), ps.sigv));
        if(s.sigu[0] > 0)
            require(isValidSignature(addr[s.u], hashChStateToSign(s), s.sigu));
        else 
            require(isValidSignature(addr[s.v], hashChStateToSign(s), s.sigv));
        chState memory os = stateForWithdl[u][v];
        if(os.v != v)
            os = symmetricChState(os);
        require(os.u == u && os.v == v);
        require(s.k > os.k);
        require(verifyStateTransition(ps, s));
        stateForWithdl[u][v] = stateForWithdl[v][u] = s;
    }

    function isValidSignature(address senderAddr, uint256 amount, bytes memory signature)
        internal
        view
        returns (bool)
    {
        bytes32 message = prefixed(keccak256(abi.encodePacked(this, amount)));
        return recoverSigner(message, signature) == senderAddr;
    }

    function splitSignature(bytes memory sig)
        internal
        pure
        returns (uint8 v, bytes32 r, bytes32 s)
    {
        require(sig.length == 65);
        assembly {
            r := mload(add(sig, 32))
            s := mload(add(sig, 64))
            v := byte(0, mload(add(sig, 96)))
        }
        return (v, r, s);
    }

    function recoverSigner(bytes32 message, bytes memory sig)
        internal
        pure
        returns (address)
    {
        (uint8 v, bytes32 r, bytes32 s) = splitSignature(sig);
        return ecrecover(message, v, r, s);
    }


    function prefixed(bytes32 hash) internal pure returns (bytes32) {
        return keccak256(abi.encodePacked("\x19Ethereum Signed Message:\n32", hash));
    }

}
