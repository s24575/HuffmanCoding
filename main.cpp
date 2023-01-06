#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

struct Node
{
    unsigned char ch = '\0';
    int freq = 0;
    Node *left = nullptr;
    Node *right = nullptr;

    ~Node()
    {
        delete left;
        delete right;
    }

    Node() = delete;
    Node(unsigned char ch, int freq) : ch(ch), freq(freq) {}
    Node(unsigned char ch, int freq, Node *left, Node *right) : ch(ch), freq(freq), left(left), right(right) {}
};

namespace pq
{
    template <class _Ty, class _Pr = std::less<_Ty>>
    class priority_queue
    {
    public:
        priority_queue() = default;

        explicit priority_queue(const _Pr &_Pred) : c(), comp(_Pred) {}

        priority_queue(std::vector<_Ty> c) : c(c)
        {
            make_heap(c);
        }

        [[nodiscard]] bool empty() const
        {
            return c.empty();
        }

        [[nodiscard]] std::size_t size() const
        {
            return c.size();
        }

        [[nodiscard]] _Ty top() const
        {
            return c.front();
        }

        void push(const _Ty &a)
        {
            c.push_back(a);
            heapifyUp();
        }

        template <class... Args>
        void emplace(Args &&...args)
        {
            c.emplace_back(std::forward<Args>(args)...);
            heapifyUp();
        }

        void pop()
        {
            std::swap(c[0], c.back());
            c.pop_back();
            heapifyDown(0);
        }

        void print()
        {
            std::cout << "Queue:" << '\n';
            const auto s = size();
            for (int i = 0; i < s; ++i)
            {
                std::cout << c[i]->ch << ' ' << c[i]->freq << '\n';
            }
            std::cout << '\n';
        }

    protected:
        [[nodiscard]] size_t getLeftChildIndex(const int index) const { return index * 2 + 1; }
        [[nodiscard]] size_t getRightChildIndex(const int index) const { return index * 2 + 2; }
        [[nodiscard]] size_t getParentIndex(const int index) const { return (index - 1) / 2; }
        [[nodiscard]] bool hasLeftChild(const size_t index) const { return getLeftChildIndex(index) < size(); }
        [[nodiscard]] bool hasRightChild(const size_t index) const { return getRightChildIndex(index) < size(); }
        [[nodiscard]] bool hasParent(const size_t index) const { return getParentIndex(index) >= 0; }
        [[nodiscard]] _Ty getLeftChild(const int index) const { return c[getLeftChildIndex(index)]; }
        [[nodiscard]] _Ty getRightChild(const int index) const { return c[getRightChildIndex(index)]; }
        [[nodiscard]] _Ty getParent(const int index) const { return c[getParentIndex(index)]; }

        void make_heap()
        {
            for (int i = size() / 2; i >= 0; i--)
            {
                heapifyDown(i);
            }
        }

        void heapifyDown(int index)
        {
            while (hasLeftChild(index))
            {
                int smallerChildIndex = getLeftChildIndex(index);
                if (hasRightChild(index) && comp(getLeftChild(index), getRightChild(index)))
                {
                    smallerChildIndex = getRightChildIndex(index);
                }

                if (comp(c[smallerChildIndex], c[index]))
                {
                    break;
                }

                std::swap(c[smallerChildIndex], c[index]);
                index = smallerChildIndex;
            }
        }

        void heapifyUp()
        {
            int index = size() - 1;
            while (index != 0 && comp(getParent(index), c[index]))
            {
                std::swap(c[index], c[getParentIndex(index)]);
                index = getParentIndex(index);
            }
        }

        std::vector<_Ty> c{};
        _Pr comp{};
    };
}

// <------------------------------------START HERE------------------------------------>

void assignCodeToHuffmanTree(Node *node, std::unordered_map<unsigned char, std::vector<bool>> &m, std::vector<bool> &code)
{
    if (!node->left && !node->right)
    {
        m[node->ch] = code;
        return;
    }

    if (node->left)
    {
        auto left_copy = code;
        left_copy.push_back(1);
        assignCodeToHuffmanTree(node->left, m, left_copy);
    }

    if (node->right)
    {
        auto right_copy = code;
        right_copy.push_back(0);
        assignCodeToHuffmanTree(node->right, m, right_copy);
    }
}

std::unordered_map<unsigned char, std::vector<bool>> generateCodes(const std::string &input)
{
    std::unordered_map<unsigned char, int> freq;
    for (const auto &ch : input)
    {
        ++freq[ch];
    }

    auto compare = [](const Node *const a, const Node *const b)
    {
        return a->freq > b->freq;
    };

    pq::priority_queue<Node *, decltype(compare)> q(compare);

    for (const auto &it : freq)
    {
        q.emplace(new Node(it.first, it.second));
    }

    while (q.size() > 1)
    {
        auto x = q.top();
        q.pop();
        auto y = q.top();
        q.pop();

        q.emplace(new Node('\0', x->freq + y->freq, x, y));
    }

    Node *root = q.top();
    std::unordered_map<unsigned char, std::vector<bool>> codes;
    std::vector<bool> code(1, 0);
    assignCodeToHuffmanTree(root, codes, code);
    delete root;

    return codes;
}

std::string convertBitsToString(const std::vector<bool> &bits)
{
    constexpr int bits_in_char = sizeof(unsigned char) * 8;
    size_t s = (bits.size() + bits_in_char - 1) / bits_in_char;
    std::string result(s, '\0');
    auto char_ptr = result.begin();
    int shift = 7;
    for (const bool bit : bits)
    {
        *char_ptr |= bit << shift;

        if (--shift < 0)
        {
            ++char_ptr;
            shift = 7;
        }
    }
    return result;
}

std::vector<bool> convertStringToBits(const std::string &str)
{
    std::vector<bool> result;
    for (const auto &ch : str)
    {
        for (int shift = 7; shift >= 0; --shift)
        {
            result.push_back(ch & 1 << shift);
        }
    }

    return result;
}

std::vector<bool> encodeString(std::string str, const std::unordered_map<unsigned char, std::vector<bool>> &codes)
{
    std::vector<bool> all_bits;

    for (const auto &ch : str)
    {
        const auto &bits = codes.at(ch);
        for (const bool bit : bits)
        {
            all_bits.push_back(bit);
        }
    }

    return all_bits;
}

std::string decodeString(const std::vector<bool> &bits, const std::unordered_map<std::vector<bool>, unsigned char> &codes)
{
    std::string decoded = "";
    std::vector<bool> curr;

    for (const bool bit : bits)
    {
        curr.push_back(bit);

        if (codes.find(curr) != codes.end())
        {
            decoded += codes.at(curr);
            curr.clear();
        }
    }

    return decoded;
}

int main()
{
    // ENCODING PART

    // open input file
    std::ifstream input("unprocessed.txt");
    if (input.fail())
    {
        std::cout << "Input file failed to open.\n";
        return -1;
    }

    // read from input file
    std::stringstream buffer;
    buffer << input.rdbuf();
    std::string unprocessed_text = buffer.str();
    if (unprocessed_text.empty())
    {
        std::cout << "The file is empty.\n";
        return 0;
    }
    // std::cout << "Input: [" << unprocessed_text << "]\n";
    input.close();

    // open output file
    std::ofstream output("processed.txt", std::ios::binary);
    if (output.fail())
    {
        std::cout << "Output file failed to open.\n";
        return -1;
    }

    // huffman coding
    auto codes = generateCodes(unprocessed_text);

    // write codes to output file
    output << codes.size() << std::endl;
    for (const auto &kv : codes)
    {
        std::string bits = "";
        for (const bool bit : kv.second)
        {
            bits += std::to_string(bit);
        }
        std::cout << bits << " " << kv.first << '\n';
        output << bits << " " << (int)kv.first << std::endl;
    }

    // write bit count and encoded string to output file
    const std::vector<bool> all_bits = encodeString(unprocessed_text, codes);
    output << all_bits.size() << std::endl;
    output << convertBitsToString(all_bits);
    output.close();

    /*for (const bool bit : all_bits) {
        std::cout << bit;
    }
    std::cout << '\n';*/

    // DECODING PART

    // read codes
    std::ifstream file("processed.txt", std::ios::in | std::ios::binary);

    int n;
    file >> n;
    std::unordered_map<std::vector<bool>, unsigned char> codes2;
    std::string key;
    int value;
    for (int i = 0; i < n; ++i)
    {
        file >> key >> value;
        // std::cout << "k:[" << key << "]v:[" << (unsigned char)value << "]\n";
        std::vector<bool> bits;
        for (const auto &ch : key)
        {
            bits.push_back((bool)(ch - '0'));
        }
        codes2[bits] = (unsigned char)value;
    }

    int bit_count;
    file >> bit_count;

    // ignore
    file.ignore(1);

    // read bits
    char c;
    std::vector<bool> bits;
    while (file.get(c))
    {
        for (int i = 7; i >= 0; i--)
        {
            bits.push_back((c >> i) & 1);
        }
    }

    // remove last few extra bits
    bits.resize(bit_count);

    /*for (const bool bit : bits) {
        std::cout << bit;
    }
    std::cout << std::endl;*/

    // check if it works
    std::string processed_text = decodeString(bits, codes2);
    std::cout << "Decoded final: [" << processed_text << "]\n";
}