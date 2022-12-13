#include "FuncMaker.hpp"

const std::vector<std::string> FunctionalTree::operations = {"+", "-", "*", "/", "%", "^", "sqrt", "sin", "cos", "tan", "cot", "asin", "acos", "atan", "acot", "log", "ln", "exp", "abs"};

FunctionalTreeNode::FunctionalTreeNode (NodeType type) : type(type), priority(0) {}
FunctionalTreeNode::~FunctionalTreeNode () {}

OperationNode::OperationNode (Operation op) : FunctionalTreeNode(NodeType::OPERATION), op(op) {}
OperationNode::~OperationNode () {}

ValueNode::ValueNode (double val) : FunctionalTreeNode(NodeType::VALUE), val(val) {}
ValueNode::~ValueNode () {}

VariableNode::VariableNode (uint64_t idx) : FunctionalTreeNode(NodeType::VARIABLE), idx(idx) {}
VariableNode::~VariableNode () {}

void FunctionalTree::inputCheck (const std::vector<std::string> &vars) const {
    if (vars.size() > VARIABLE_LIMIT) {
        throw std::logic_error("Operation \"inputCheck\": count of vars limited by " + std::to_string(VARIABLE_LIMIT));
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        std::vector<std::string>::const_iterator it = std::find(operations.cbegin(), operations.cend(), vars[i]);
        uint64_t idx = std::distance(operations.cbegin(), it);
        if (idx != operations.size()) {
            throw std::logic_error("Operation \"inputCheck\": used var name \"" + vars[i] + "\" which is a name of operation");
        }
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        uint64_t count = std::count(vars.cbegin(), vars.cend(), vars[i]);
        if (count != 1) {
            throw std::logic_error("Operation \"inputCheck\": var name \"" + vars[i] + "\" used more than 1 time");
        }
    }
}

std::string FunctionalTree::readOperation (const std::string &func, uint64_t &i) const {
    std::string str;
    bool flag = true;
    while (i < func.size() && flag) {
        switch (func[i]) {
            case '+':
            case '-':
            case '*':
            case '/':
            case '%':
            case '^':
                str += func[i];
                ++i;
                break;
            default:
                flag = false;
                break;
        }
    }
    return str;
}

std::string FunctionalTree::readWord (const std::string &func, uint64_t &i) const {
    std::string str;
    //words
    while (((func[i] >= 'A' && func[i] <= 'Z') || (func[i] >= 'a' && func[i] <= 'z') || func[i] == '\'' || func[i] == '\"') && i < func.size()) {
        str += func[i];
        ++i;
    }
    //+,-,*,/,%,^
    if (str.empty() && i < func.size()) {
        str = readOperation(func, i);
    }
    return str;
}

double FunctionalTree::readNumber (const std::string &func, uint64_t &i) const {
    std::string str;
    while (((func[i] >= '0' && func[i] <= '9') || func[i] == '.') && i < func.size()) {
        str += func[i];
        ++i;
    }
    return std::atof(str.c_str());
}

std::string FunctionalTree::readInbrace (const std::string &func, uint64_t &i) const {
    uint64_t braceCount = 1;
    std::string str;
    ++i;
    while (i < func.size()) {
        if (func[i] == '(') {
            ++braceCount;
        } else if (func[i] == ')') {
            --braceCount;
        }
        if (braceCount != 0) {
            str += func[i];
        } else {
            break;
        }
        ++i;
    }
    ++i;
    if (braceCount != 0) {
        throw std::out_of_range("Operation \"readInbrace\": out of range. Incorrect placement of brackets");
    }
    return str;
}

Operation FunctionalTree::getOperation (const std::string &str) const {
    Operation op = Operation::NOT_AN_OPERATION;
    if (str == "**") {
        return Operation::POW;
    }
    for (uint64_t i = 0; i < operations.size(); ++i) {
        if (str == operations[i]) {
            op = static_cast<Operation>(i);
            break;
        }
    }
    return op;
}

uint64_t FunctionalTree::getPriority (Operation op) const {
    switch (op) {
        case Operation::PLUS:
        case Operation::MINUS:
            return 4;
        case Operation::MUL:
        case Operation::DIV:
        case Operation::MOD:
            return 3;
        case Operation::POW:
            return 2;
        default:
            return 1;
    }
}

double FunctionalTree::useOperation (Operation op, double x, double y) const {
    switch (op) {
        case Operation::PLUS:
            return x + y;
        case Operation::MINUS:
            return x - y;
        case Operation::MUL:
            return x * y;
        case Operation::DIV:
            return x / y;
        case Operation::MOD:
            return std::fmod(x, y);
        case Operation::POW:
            return std::pow(x, y);
        case Operation::SQRT:
            return std::sqrt(y);
        case Operation::SIN:
            return std::sin(y);
        case Operation::COS:
            return std::cos(y);
        case Operation::TAN:
            return std::tan(y);
        case Operation::CTG:
            return 1.0 / std::tan(y);
        case Operation::ASIN:
            return std::asin(y);
        case Operation::ACOS:
            return std::acos(y);
        case Operation::ATAN:
            return std::atan(y);
        case Operation::ACOT:
            return std::acos(-1) / 2 -  std::atan(y);
        case Operation::LOG:
            return std::log10(y);
        case Operation::LN:
            return std::log(y);
        case Operation::EXP:
            return std::exp(y);
        case Operation::ABS:
            return std::abs(y);
        default:
            return 0;
    }
}

double FunctionalTree::getVal (const NodePtr &node, const std::vector<double> &X) const {
    if (!node->left.get() && !node->right.get()) {
        if (node->type == NodeType::VALUE) {
            ValueNode* value = (ValueNode*) node.get();
            return value->val;
        } else {
            VariableNode* var = (VariableNode*) node.get();
            return X[var->idx];
        }
    }
    OperationNode* operation = (OperationNode*) node.get();
    if (node->left && !node->right) {
        //return useOperation(operation->op, getVal(operation->left, X), 0);
        return useOperation(operation->op, 0, getVal(operation->left, X));
    }
    if (!node->left && node->right) {
        //return useOperation(operation->op, getVal(operation->right, X), 0);
        return useOperation(operation->op, 0, getVal(operation->right, X));
    }
    double a = getVal(node->left, X);
    double b = getVal(node->right, X);
    return useOperation(operation->op, a, b);
}

void FunctionalTree::addToTree (NodePtr &tree, NodePtr &toAdd) {
    if (!tree.get()) {
        tree.swap(toAdd);
        return;
    }
    //switch (toAdd->type) {
    //    case NodeType::OPERATION:
            if (tree->priority > toAdd->priority) {
                if (!tree->left) {
                    tree->left.swap(toAdd);
                } else if (!tree->right) {
                    tree->right.swap(toAdd);
                } else if (tree->right->priority > toAdd->priority) {
                    // if (!tree->right->left) {
                    //     tree->right->left.swap(toAdd);
                    // } else {
                    //     tree->right->right.swap(toAdd);
                    // }
                    addToTree(tree->right, toAdd);
                } else {
                    toAdd->left.swap(tree->right);
                    tree->right.swap(toAdd);
                }
            } else {
                toAdd->left.swap(tree);
                tree.swap(toAdd);
            }
    //         break;
    //     case NodeType::VALUE:
    //     case NodeType::VARIABLE:
    //         if (!tree->left.get()) {
    //             tree->left.swap(toAdd);
    //         } else {
    //             auto tmp = tree.get();
    //             while (tmp->right.get()) {
    //                 tmp = tmp->right.get();
    //             }
    //             if (!tmp->left.get()) {
    //                 tmp->left.swap(toAdd);
    //             } else {
    //                 tmp->right.swap(toAdd);
    //             }
    //         }
    //         break;
    //     default:
    //         break;
    // }
}

FunctionalTree::NodePtr FunctionalTree::buildTree (const std::string &func, const std::vector<std::string> &vars) {
    //std::cout << "Our func: " << func << "\n";
    std::string tmp;
    uint64_t i = 0;

    double num;
    Operation op;
    uint64_t idx;
    NodePtr currentNode, node;
    while (i < func.size()) {
        //std::cout << "new step: i = " << i << "\n";
        if (func[i] == ' ') {
            ++i;
            continue;
        }
        if ((func[i] >= '0' && func[i] <= '9') || func[i] == '.') {
            num = readNumber(func, i);
            currentNode.reset(new ValueNode(num));
            //std::cout << "got num: " << number << "!\n";
        } else if (func[i] == '(') {
            tmp = readInbrace(func, i);
            //std::cout << "inbraced: " << tmp << "! Building new tree\n";
            currentNode = buildTree(tmp, vars);
            //if (currentNode->type == NodeType::OPERATION) {
                //((OperationNode *)currentNode.get())->priority = 0;
            currentNode->priority = 0;
           //}
        } else {
            //tmp = readOperation(func, i);
            tmp = readWord(func, i);
            op = getOperation(tmp);
            if (op == Operation::NOT_AN_OPERATION) {
                std::vector<std::string>::const_iterator it = std::find(vars.cbegin(), vars.cend(), tmp);
                idx = std::distance(vars.begin(), it);
                if (idx == vars.size()) {
                    throw std::logic_error("Operation \"buildTree\": var \"" + tmp + "\" not found in var list");;
                }
                currentNode.reset(new VariableNode(idx));
                //std::cout << "got var: " << tmp << "!\n";
            } else {
                currentNode.reset(new OperationNode(op));
                currentNode->priority = getPriority(op);
                //std::cout << "got operation: " << tmp << "!\n";
                //((OperationNode *)currentNode.get())->priority = getPriority(op);
            }
        }
        //std::cout << "read: " << tmp << "!\n";
        addToTree(node, currentNode);
    }
    return node;
}

void FunctionalTree::printTree (const NodePtr &node) const {
    if (!node) {
        return;
    }
    switch (node->type) {
        case NodeType::VALUE:
            std::cout << ((ValueNode *) node.get())->val;
            break;
        case NodeType::VARIABLE:
            std::cout << "x" << ((VariableNode *) node.get())->idx + 1;;
            break;
        case NodeType::OPERATION:
            std::cout << operations[static_cast<uint64_t>(((OperationNode *) node.get())->op)];
            break;
        default:
            break;
    }
    std::cout << "\n";
    if (node->left) {
        printTree(node->left);
    } else {
        std::cout << "no left node\n";
    }
    if (node->right) {
        printTree(node->right);
    } else {
        std::cout << "no right node\n";
    }
}

void FunctionalTree::printFunc (const NodePtr &node) const {
    if (!node) {
        return;
    }
    if (!node->left && !node->right) {
        switch (node->type) {
            case NodeType::VALUE:
                std::cout << ((ValueNode *) node.get())->val;
                break;
            case NodeType::VARIABLE:
                std::cout << "x" << ((VariableNode *) node.get())->idx + 1;
                break;
            default:
                break;
        }
    } else if (!node->left != !node->right) {
        std::cout << operations[static_cast<uint64_t>(((OperationNode *) node.get())->op)];
        if (((OperationNode *) node.get())->op == Operation::MINUS) {
            printFunc(node->left ? node->left : node->right);
        } else {
            std::cout << "(";
            printFunc(node->left ? node->left : node->right);
            std::cout << ")";
        }
    } else {
        if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
            std::cout << "(";
            printFunc(node->left);
            std::cout << ")";
        } else {
            printFunc(node->left);
        }
        std::cout << operations[static_cast<uint64_t>(((OperationNode *) node.get())->op)];
        if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
            std::cout << "(";
            printFunc(node->right);
            std::cout << ")";
        } else {
            printFunc(node->right);
        }
    }
}

FunctionalTree::NodePtr FunctionalTree::copyTree (const NodePtr &toCopy) const {
    if (!toCopy) {
        return nullptr;
    }
    //NodePtr node(new FunctionalTreeNode(toCopy->type));
    NodePtr node;
    auto tmp1 = (const OperationNode *) toCopy.get();
    auto tmp2 = (const ValueNode *) toCopy.get();
    auto tmp3 = (const VariableNode *) toCopy.get();
    //FunctionalTreeNode *tmp;
    switch (toCopy->type) {
        case NodeType::OPERATION:
            //auto tmp1 = (const OperationNode *) toCopy.get();
            node.reset(new OperationNode(tmp1->op));
            break;
        case NodeType::VALUE:
            //auto tmp2 = (const ValueNode *) toCopy.get();
            node.reset(new ValueNode(tmp2->val));
            break;
        case NodeType::VARIABLE:
            //auto tmp3 = (const VariableNode *) toCopy.get();
            node.reset(new VariableNode(tmp3->idx));
            break;
        default:
            break;
    }
    node->priority = toCopy->priority;
    node->left = copyTree(toCopy->left);
    node->right = copyTree(toCopy->right);
    return node;
}

FunctionalTree::FunctionalTree (const NodePtr &node) {
    root = copyTree(node);
}

FunctionalTree::FunctionalTree (NodePtr &&tree) {
    root = std::move(tree);
}

FunctionalTree::FunctionalTree () {}

FunctionalTree::FunctionalTree (const std::string &func, const std::vector<std::string> &vars) {
    reset(func, vars);
}

// FunctionalTree::FunctionalTree (const std::string &func, const std::string &var) {
//     FunctionalTree(func, {var});
// }

FunctionalTree::FunctionalTree (const FunctionalTree &tree) {
    root = copyTree(tree.root);
    //vars = tree.vars;
}

FunctionalTree::FunctionalTree (FunctionalTree &&tree) {
    root = std::move(tree.root);
    //vars = std::move(tree.vars);
}

FunctionalTree::~FunctionalTree () {}

void FunctionalTree::reset (const std::string &func, const std::vector<std::string> &vars) {
    inputCheck(vars);
    //this->vars = vars;
    root = buildTree(func, vars);
}

double FunctionalTree::func (double x) const {
    return getVal(root, {x});
}

double FunctionalTree::func (const std::vector<double> &X) const {
    return getVal(root, X);
}

FunctionalTree FunctionalTree::getCoeff (uint64_t idx) const {
    auto tmp = root.get();
    while (tmp->left && tmp->type == NodeType::OPERATION) {
        if (tmp->right) {
            if (tmp->right->type == NodeType::VARIABLE) {
                auto n = (VariableNode *) tmp->right.get();
                auto o = (OperationNode *) tmp;
                //std::cout << "here\n";
                if (n->idx == idx) {
                    if (o->op == Operation::MINUS) {
                        return FunctionalTree("-1", {});
                    }
                    if (o->op == Operation::PLUS) {
                        return FunctionalTree("1", {});
                    }
                    if (o->op == Operation::MUL) {
                        return FunctionalTree(tmp->left);
                    }
                }
            }
            if (tmp->right->right && tmp->right->right->type == NodeType::VARIABLE) {
                auto n = (VariableNode *) tmp->right->right.get();
                auto o = (OperationNode *) tmp;
                if (n->idx == idx) {
                    //std::cout << "got it!\n";
                    FunctionalTree tr(tmp->right->left);
                    if (o->op == Operation::MINUS) {
                        //std::cout << "got it! again...\n";
                        tr.root->priority = 0;
                        auto node = tr.buildTree("-", {});
                        tr.addToTree(tr.root, node);
                    }
                    return tr;
                }
            }
        }
        tmp = tmp->left.get();
    }
    if (tmp->type == NodeType::VARIABLE) {
        auto n = (VariableNode *) tmp;
        //std::cout << "i: " << idx << " " << n->idx << "\n";
        if (n->idx == idx) {
            return FunctionalTree("1", {});
        }
    }
    return FunctionalTree("0", {});
    // if (tmp->type != NodeType::OPERATION || !tmp->right || !tmp->right->left) {
    //     return FunctionalTree("0", {});
    // }
    // if (tmp->right->type == NodeType::VARIABLE) {
    //     auto n = (VariableNode *) tmp->right.get();
    // }
    // return FunctionalTree(tmp->right->left);
}

FunctionalTree FunctionalTree::getDiv () const {
    if (root->type == NodeType::OPERATION) {
        return FunctionalTree(root->right);
    } else {
        return FunctionalTree();
    }
}

void FunctionalTree::printTree () const {
    printTree(root);
}

void FunctionalTree::printFunc () const {
    printFunc(root);
}

//void FunctionalTree::simplify () {}

FunctionalTree &FunctionalTree::operator= (const FunctionalTree &tree) {
    if (this == &tree) {
        return *this;
    }
    root = copyTree(tree.root);
    //vars = tree.vars;
    return *this;
}

double FunctionalTree::operator() (double x) const {
    return getVal(root, {x});
}
double FunctionalTree::operator() (const std::vector<double> &X) const {
    return getVal(root, X);
}

FunctionalTree &FunctionalTree::operator= (FunctionalTree &&tree) {
    if (this == &tree) {
        return *this;
    }
    root = std::move(tree.root);
    //vars = std::move(tree.vars);
    return *this;
}