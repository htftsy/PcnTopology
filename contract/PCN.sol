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
    uint32 public channelPerIndex;
    uint32 public custimizedChannelPerIndex;
    uint32 public dutyChannelPerIndex;
    uint256 public amountPerChannel;
    uint32[] public remainedCustimizedChannels;

    struct matrice {
        uint32 a;
        uint32 b;
        uint32 c;
        uint32 d;
    }

    matrice[] public PGLMembers;
    matrice[] public SMembers;
    mapping(uint32 => uint32) public matriceToIndex; 
    mapping(uint32 => uint32) public iterNumToIndex;
    mapping(uint32 => mapping(uint32 => uint32)) public adjacentIndex;
    mapping(uint32 => mapping(uint32 => uint256)) public capacity;

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

    function determineIterOrder() public {
        iterNumToIndex[0] = 0;
        uint32 front = 0;
        uint32 rear = 0;
        uint32[] memory qu;
        qu = new uint32[] (N + 1);
        qu[0] = 0;
        bool[] memory vis;
        vis = new bool[] (N + 1);
        vis[0] = true;
        for(uint32 i = 1; i < N; i++)
            vis[i] = false;
        uint32 target;
        for( ; front <= rear; front ++) {
            uint32 f = qu[front];
            for(uint32 j = 0; j <= P; j++) {
                target = adjacentIndex[f][j];
                if(! vis[target]) {
                    qu[++ rear] = target;
                    vis[target] = true;
                    iterNumToIndex[rear] = target;
                }
            }
        }
    }

    function setupLPS(uint32 _P, uint32 _Q) public {
        P = _P;
        Q = _Q;
        
        genPGLMembers();
        genSMembers();

        matrice memory t;

        for(uint32 i = 0; i < N; i++)
            matriceToIndex[serializeMatriceToInt(PGLMembers[i])] = i;

        for(uint32 i = 0; i < N; i++) {
            for(uint32 j = 0; j <= uint32(P); j++) {
                t = matriceMul(PGLMembers[i], SMembers[j]);
                uint32 k = matriceToIndex[serializeMatriceToInt(t)];
                adjacentIndex[i][j] = k;
            }
        }

        determineIterOrder();
    }

    constructor (
        uint32 _P, 
        uint32 _Q, 
        uint32 custimized_channel_per_index, // this is (1 - gamma) * m in the paper
        uint256 amount_per_channel
    )  
    {
        uint32 channel_per_index = _P + 1;
        require(custimized_channel_per_index <= channel_per_index);

        N = (_Q ** 3) - _Q;    // `**' means power for Solidity

        iter = 0;
        adj = new uint32[][](N);  
        addr = new address payable[](N);

        amountPerChannel = amount_per_channel;
        channelPerIndex = channel_per_index;
        custimizedChannelPerIndex = custimized_channel_per_index;
        dutyChannelPerIndex = channelPerIndex - custimizedChannelPerIndex;

        remainedCustimizedChannels = new uint32[](N);

        setupLPS(_P, _Q);
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
        require(msg.value >= amountPerChannel * uint256(channelPerIndex));   // deposits for both duty channels and typical channels

        uint32 index = iterNumToIndex[iter];
        iter ++;
        require(corrdinateToIndex(indexToCoordinate(index)) == index);

        addr[index] = payable(msg.sender);

        adj[index] = new uint32[](dutyChannelPerIndex);

        for(uint32 k = 0; k < dutyChannelPerIndex; k++)
                adj[index][k] = adjacentIndex[index][k];

        for(uint32 i = 0; i < dutyChannelPerIndex; i++)
            capacity[index][uint32(adj[index][uint32(i)])] = amountPerChannel;

        remainedCustimizedChannels[index] = custimizedChannelPerIndex;

        return (index, adj[index]);
    }

    function addCustimizedChannel (uint32 index, uint32 indexTo, bytes memory signature) 
        public 
    {
        require(addr[index] == msg.sender);
        require(remainedCustimizedChannels[index] > 0);
        require(isValidSignature(addr[indexTo], uint256(index), signature));
        adj[index].push(indexTo);
        capacity[index][indexTo] = amountPerChannel;
        remainedCustimizedChannels[index] --;
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

    function queryPrescribedInwardAdjacentIndices (uint32 index) 
        public 
        view
        returns (uint32[] memory res) 
    {
        res = new uint32[](dutyChannelPerIndex);
        uint32 count = 0;

        for(uint32 k = 0; k < dutyChannelPerIndex; k++) {
            res[count] = adjacentIndex[index][k];
            count ++;
        }

        return res;
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

    function liquidate (uint32 _indexSender, uint32 _indexReceiver, address senderAddr, uint256 amount, bytes memory signature) 
        public 
    {
        uint32 indexSender = _indexSender;
        uint32 indexReceiver = _indexReceiver;
        require(amount <= capacity[indexSender][indexReceiver]);
        require(addr[indexSender] == senderAddr);
        require(addr[indexReceiver] == msg.sender);
        require(isValidSignature(senderAddr, amount, signature));
        uint32 len = uint32(adj[indexSender].length);
        bool found = false;
        for(uint i = 0; i < len; i++) {
            if(uint32(adj[indexSender][i]) == indexReceiver) {
                found = true;
                adj[indexSender][i] = adj[indexSender][len - 1];
                adj[indexSender].pop(); // channel removed thereafter
            }
        }
        require(found); 
        addr[indexReceiver].transfer(amount);
        capacity[indexSender][indexReceiver] -= amount;
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