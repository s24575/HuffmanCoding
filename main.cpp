#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

struct HuffmanNode {
    char ch = '\0';
    int freq = 0;
    HuffmanNode* parent = nullptr;
    HuffmanNode* left = nullptr;
    HuffmanNode* right = nullptr;

    ~HuffmanNode() {
        delete left;
        delete right;
    }

    HuffmanNode() = delete;
    HuffmanNode(int freq) : freq(freq) {}
    HuffmanNode(char ch, int freq) : ch(ch), freq(freq) {}

    constexpr bool operator<(const HuffmanNode& _Right) const noexcept {
        return this->freq < _Right.freq;
    }
};

namespace pq {

    template<typename _Ty = void>
    struct comp {
        constexpr bool operator()(const _Ty& _Left, const _Ty& _Right) const noexcept {
            return _Left < _Right;
        }
    };

    template<class _Ty, class _Pr = comp<_Ty>>
    class priority_queue {
    public:
        priority_queue() = default;

        priority_queue(std::vector<_Ty> c) : c(c) {
            make_heap(c);
        }

        [[nodiscard]] bool empty() const {
            return c.empty();
        }

        [[nodiscard]] std::size_t size() const {
            return c.size();
        }

        [[nodiscard]] _Ty top() const {
            return c.front();
        }

        void push(const _Ty& a) {
            c.push_back(a);
            heapifyUp();
        }

        template <class... Args>
        void emplace(Args&&... args) {
            c.emplace_back(std::forward<Args>(args)...);
            heapifyUp();
        }

        void pop() {
            std::swap(c[0], c.back());
            c.pop_back();
            heapifyDown(0);
        }

    protected:
        [[nodiscard]] bool hasLeftChild(const size_t index) const { return getLeftChildIndex(index) < size(); }
        [[nodiscard]] bool hasRightChild(const size_t index) const { return getRightChildIndex(index) < size(); }
        [[nodiscard]] bool hasParent(const size_t index) const { return getParentIndex(index) >= 0; }
        [[nodiscard]] size_t getLeftChildIndex(const int index) const { return index * 2 + 1; }
        [[nodiscard]] size_t getRightChildIndex(const int index) const { return index * 2 + 2; }
        [[nodiscard]] size_t getParentIndex(const int index) const { return (index - 1) / 2; }
        [[nodiscard]] _Ty getLeftChild(const int index) const { return c[getLeftChildIndex(index)]; }
        [[nodiscard]] _Ty getRightChild(const int index) const { return c[getRightChildIndex(index)]; }
        [[nodiscard]] _Ty getParent(const int index) const { return c[getParentIndex(index)]; }

        void make_heap() {
            for (int i = size() / 2; i >= 0; i--) {
                heapifyDown(i);
            }
        }

        void heapifyDown(int index) {
            while (hasLeftChild(index)) {
                int smallerChildIndex = getLeftChildIndex(index);
                if (hasRightChild(index) && comp(getRightChild(index), getLeftChild(index))) {
                    smallerChildIndex = getRightChildIndex(index);
                }

                if (comp(c[index], c[smallerChildIndex])){
                    break;
                }
                else {
                    std::swap(c[index], c[smallerChildIndex]);
                    index = smallerChildIndex;
                }
            }
        }

        void heapifyUp() {
            int index = size() - 1;
            while (hasParent(index) && comp(getParent(index), c[index])) {
                int parentIndex = getParentIndex(index);
                std::swap(c[index], c[parentIndex]);
                index = parentIndex;
            }
        }

        std::vector<_Ty> c{};
        _Pr comp{};
    };
}

void assignCodeToHuffmanTree(HuffmanNode* node, std::unordered_map<char, std::string>& m, std::string code = "") {
    if (node->ch) {
        m[node->ch] = code;
    }
    else {
        if (node->left) {
            assignCodeToHuffmanTree(node->left, m, code + "1");
        }
        if (node->right) {
            assignCodeToHuffmanTree(node->right, m, code + "0");
        }
    }
}

std::unordered_map<char, std::string> generateCodes(const std::string& input) {
    std::unordered_map<char, int> freq;
    for (const auto& i : input) {
        ++freq[i];
    }

    pq::priority_queue<HuffmanNode*> q;

    for (const auto& it : freq) {
        q.emplace(new HuffmanNode(it.first, it.second));
    }

    while (q.size() > 1) {
        auto x = q.top();
        q.pop();
        auto y = q.top();
        q.pop();
        auto parent = new HuffmanNode(x->freq + y->freq);
        parent->left = x;
        parent->right = y;
        x->parent = parent;
        y->parent = parent;
        q.push(parent);
    }

    HuffmanNode* root = q.top();
    std::unordered_map<char, std::string> codes;
    assignCodeToHuffmanTree(root, codes);
    delete root;

    return codes;
}

std::string encodeString(std::string str, const std::unordered_map<char, std::string>& codes) {
    std::string result = "";
    for (const auto& ch : str) {
        result += codes.at(ch);
    }
    return result;
}

std::string decodeString(std::string file_name) {
    std::ifstream file(file_name);
    int n;
    file >> n;
    std::unordered_map<std::string, char> codes;
    std::string key;
    char value;
    for (int i = 0; i < n; ++i) {
        file >> key >> value;
        codes[key] = value;
    }

    std::string text, encoded;
    file >> text;
    std::cout << text << '\n';
    std::string curr = "";
    for (const auto& ch : text) {
        curr += ch;
        auto elem = codes.find(curr);
        if (elem != codes.end()) {
            encoded += codes[curr];
            curr = "";
        }
    }
    return encoded;
}

int main() {
    std::ifstream file("unprocessed.txt");
    if (file.fail()) {
        std::cout << "File failed to open.\n";
        return -1;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string unprocessed_text = buffer.str();
    std::cout << '[' << unprocessed_text << ']' << '\n';
    file.close();

    auto codes = generateCodes(unprocessed_text);
    for (const auto& kv : codes) {
        std::cout << "[" << kv.first << "]" << ' ' << kv.second << '\n';
    }

    std::ofstream out("processed.txt");
    out << codes.size() << '\n';
    for (const auto& kv : codes) {
        out << kv.second << " " << kv.first << std::endl;
    }
    std::string processed_text = encodeString(unprocessed_text, codes);
    out << processed_text;
    out.close();

    std::string reverse = decodeString("processed.txt");
    std::cout << '[' << reverse << ']' << '\n';
}