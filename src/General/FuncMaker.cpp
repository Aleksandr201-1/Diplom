#include "FuncMaker.hpp"

//var style

OperationStruct::OperationStruct (const std::vector<std::string> &op_str, const std::function<double (double, double)> &func, uint64_t priority) 
: op_str(op_str), func(func), priority(priority) {}

OperationStruct::~OperationStruct () {}

ConstantValue::ConstantValue (const std::string &val_name, double val) : val_name(val_name), val(val) {}
ConstantValue::~ConstantValue () {}

const std::vector<OperationStruct> FunctionalTree::operations = {
    {{"+"},                        [] (double x, double y) {return x + y;},                                4},
    {{"-"},                        [] (double x, double y) {return x - y;},                                4},
    {{"*"},                        [] (double x, double y) {return x * y;},                                3},
    {{"/"},                        [] (double x, double y) {return x / y;},                                3},
    {{"%"},                        [] (double x, double y) {return std::fmod(x, y);},                      3},
    {{"^", "**"},                  [] (double x, double y) {return std::pow(x, y);},                       2},
    {{"sqrt"},                     [] (double x, double y) {return std::sqrt(x);},                         1},
    {{"sin"},                      [] (double x, double y) {return std::sin(x);},                          1},
    //{{"QWERTY"},                   [] (double x, double y) {return std::sin(x) + std::cos(y);},            1},
    {{"cos"},                      [] (double x, double y) {return std::cos(x);},                          1},
    {{"tan", "tg"},                [] (double x, double y) {return std::tan(x);},                          1},
    {{"cot", "ctg"},               [] (double x, double y) {return 1.0 / std::tan(x);},                    1},
    {{"sinh"},                     [] (double x, double y) {return std::sinh(x);},                         1},
    {{"cosh"},                     [] (double x, double y) {return std::cosh(x);},                         1},
    {{"tanh", "tgh"},              [] (double x, double y) {return std::tanh(x);},                         1},
    {{"coth", "ctgh"},             [] (double x, double y) {return 1.0 / std::tanh(x);},                   1},
    {{"arcsin", "asin"},           [] (double x, double y) {return std::asin(x);},                         1},
    {{"arccos", "acos"},           [] (double x, double y) {return std::acos(x);},                         1},
    {{"arctan", "arctg", "atan"},  [] (double x, double y) {return std::atan(x);},                         1},
    {{"arccot", "arcctg", "acot"}, [] (double x, double y) {return std::acos(-1.0) / 2.0 - std::atan(x);}, 1},
    {{"log"},                      [] (double x, double y) {return std::log10(x);},                        1},
    {{"ln"},                       [] (double x, double y) {return std::log(x);},                          1},
    {{"exp"},                      [] (double x, double y) {return std::exp(x);},                          1},
    {{"abs"},                      [] (double x, double y) {return std::abs(x);},                          1},
    {{"sign"},                     [] (double x, double y) {return x > 0.0 ? 1 : -1;},                     1}
};

const std::vector<ConstantValue> FunctionalTree::const_val = {
    {"pi", std::acos(-1.0)},
    {"e", std::exp(1.0)}
};

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
        for (uint64_t j = 0; j < operations.size(); ++j) {
            auto it = std::find(operations[j].op_str.cbegin(), operations[j].op_str.cend(), vars[i]);
            uint64_t idx = std::distance(operations[j].op_str.cbegin(), it);
            if (idx != operations[j].op_str.size()) {
                throw std::logic_error("Operation \"inputCheck\": used var name \"" + vars[i] + "\" which is a name of operation");
            }
        }
    }
    for (uint64_t i = 0; i < const_val.size(); ++i) {
        auto it = std::find(vars.cbegin(), vars.cend(), const_val[i].val_name);
        uint64_t idx = std::distance(vars.cbegin(), it);
        if (idx != vars.size()) {
            throw std::logic_error("Operation \"inputCheck\": used var name \"" + vars[idx] + "\" which is a name of constant");
        }
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        uint64_t count = std::count(vars.cbegin() + i, vars.cend(), vars[i]);
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
    for (uint64_t i = 0; i < operations.size(); ++i) {
        for (uint64_t j = 0; j < operations[i].op_str.size(); ++j) {
            if (str == operations[i].op_str[j]) {
                op = static_cast<Operation>(i);
                break;
            }
        }
    }
    return op;
}

double FunctionalTree::getConstant (const std::string &str) const {
    double constant = 0;
    for (uint64_t i = 0; i < const_val.size(); ++i) {
        if (const_val[i].val_name == str) {
            constant = const_val[i].val;
            break;
        }
    }
    return constant;
}

uint64_t FunctionalTree::getPriority (Operation op) const {
    uint64_t idx = static_cast<uint64_t>(op);
    return operations[idx].priority;
}

double FunctionalTree::useOperation (Operation op, double x, double y) const {
    uint64_t idx = static_cast<uint64_t>(op);
    return operations[idx].func(x, y);
}

double FunctionalTree::getVal (const NodePtr &node, const std::vector<double> &X) const {
    if (!node->left.get() && !node->right.get()) {
        if (node->type == NodeType::VALUE) {
            return static_cast<ValueNode *>(node.get())->val;
        } else {
            return X[static_cast<VariableNode *>(node.get())->idx];
        }
    }
    OperationNode *operation = static_cast<OperationNode *>(node.get());
    if (node->left && !node->right) {
        return useOperation(operation->op, getVal(operation->left, X), 0);
    }
    if (!node->left && node->right) {
        return useOperation(operation->op, getVal(operation->right, X), 0);
    }
    double a = getVal(node->left, X);
    double b = getVal(node->right, X);
    return useOperation(operation->op, a, b);
}

void FunctionalTree::addToTree (NodePtr &tree, NodePtr &node) {
    if (!tree.get()) {
        if (node->type == NodeType::OPERATION && node->priority != 0) {
            OperationNode *operation = static_cast<OperationNode *>(node.get());
            if (operation->op == Operation::MINUS) {
                NodePtr zero = std::make_unique<ValueNode>(0);
                node->left = std::move(zero);
            }
        }
        tree = std::move(node);
        return;
    }
    if (tree->priority > node->priority) {
        if (!tree->left) {
            tree->left = std::move(node);
        } else if (!tree->right) {
            tree->right = std::move(node);
        } else if (tree->right->priority > node->priority) {
            addToTree(tree->right, node);
        } else {
            node->left = std::move(tree->right);
            tree->right = std::move(node);
        }
    } else {
        node->left = std::move(tree);
        tree = std::move(node);
    }
}

FunctionalTree::NodePtr FunctionalTree::buildTree (const std::string &func) {
    std::string tmp;
    uint64_t i = 0;

    double num;
    Operation op;
    uint64_t idx;
    NodePtr current, node;
    while (i < func.size()) {
        if (func[i] == ' ') {
            ++i;
            continue;
        }
        if ((func[i] >= '0' && func[i] <= '9') || func[i] == '.') {
            num = readNumber(func, i);
            current = std::make_unique<ValueNode>(num);
        } else if (func[i] == '(') {
            tmp = readInbrace(func, i);
            current = buildTree(tmp);
            current->priority = 0;
        } else {
            tmp = readWord(func, i);
            op = getOperation(tmp);
            if (op == Operation::NOT_AN_OPERATION) {
                double constant = getConstant(tmp);
                if (constant == 0.0) {
                    auto it = std::find(vars.cbegin(), vars.cend(), tmp);
                    idx = std::distance(vars.cbegin(), it);
                    if (idx == vars.size()) {
                        throw std::logic_error("Operation \"buildTree\": var \"" + tmp + "\" not found in var list");;
                    }
                    current = std::make_unique<VariableNode>(idx);
                } else {
                    current = std::make_unique<ValueNode>(constant);
                    
                }
            } else {
                current = std::make_unique<OperationNode>(op);
                current->priority = getPriority(op);
            }
        }
        addToTree(node, current);
    }
    return node;
}

void FunctionalTree::printTree (const NodePtr &node, std::ostream &out) const {
    if (!node) {
        return;
    }
    Operation op;
    switch (node->type) {
        case NodeType::VALUE:
            out << static_cast<ValueNode *>(node.get())->val;
            break;
        case NodeType::VARIABLE:
            out << vars[static_cast<VariableNode *>(node.get())->idx];
            break;
        case NodeType::OPERATION:
            op = static_cast<OperationNode *>(node.get())->op;
            out << operations[static_cast<uint64_t>(op)].op_str[0];
            break;
        default:
            break;
    }
    out << "\n";
    if (node->left) {
        printTree(node->left, out);
    } else {
        out << "no left node\n";
    }
    if (node->right) {
        printTree(node->right, out);
    } else {
        out << "no right node\n";
    }
}

void FunctionalTree::printFunc (const NodePtr &node, std::ostream &out) const {
    if (!node) {
        return;
    }
    if (!node->left && !node->right) {
        switch (node->type) {
            case NodeType::VALUE:
                out << static_cast<ValueNode *>(node.get())->val;
                break;
            case NodeType::VARIABLE:
                out << vars[static_cast<VariableNode *>(node.get())->idx];
                break;
            default:
                break;
        }
    } else if (!node->left != !node->right) {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        out << operations[static_cast<uint64_t>(op)].op_str[0];
        if (op == Operation::MINUS) {
            printFunc(node->left ? node->left : node->right, out);
        } else {
            out << "(";
            printFunc(node->left ? node->left : node->right, out);
            out << ")";
        }
    } else {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
            out << "(";
            printFunc(node->left, out);
            out << ")";
        } else {
            printFunc(node->left, out);
        }
        out << operations[static_cast<uint64_t>(op)].op_str[0];
        if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
            out << "(";
            printFunc(node->right, out);
            out << ")";
        } else {
            printFunc(node->right, out);
        }
    }
}

void FunctionalTree::toStringDefault (const NodePtr &node, std::string &str) const {
    if (!node) {
        return;
    }
    std::ostringstream str_stream;
    if (!node->left && !node->right) {
        switch (node->type) {
            case NodeType::VALUE:
                str_stream << static_cast<ValueNode *>(node.get())->val;
                str += str_stream.str();
                break;
            case NodeType::VARIABLE:
                str += vars[static_cast<VariableNode *>(node.get())->idx];
                break;
            default:
                break;
        }
    } else if (!node->left != !node->right) {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        str += operations[static_cast<uint64_t>(op)].op_str[0];
        if (op == Operation::MINUS) {
            toStringDefault(node->left ? node->left : node->right, str);
        } else {
            str += "(";
            toStringDefault(node->left ? node->left : node->right, str);
            str += ")";
        }
    } else {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
            str += "(";
            toStringDefault(node->left, str);
            str += ")";
        } else {
            toStringDefault(node->left, str);
        }
        str += " " + operations[static_cast<uint64_t>(op)].op_str[0] + " ";
        if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
            str += "(";
            toStringDefault(node->right, str);
            str += ")";
        } else {
            toStringDefault(node->right, str);
        }
    }
}

//void FunctionalTree::toStringGNUPlot (const NodePtr &node, std::string &str) const {}

void FunctionalTree::toStringLatex (const NodePtr &node, std::string &str) const {
    if (!node) {
        return;
    }
    std::ostringstream str_stream;
    if (!node->left && !node->right) {
        switch (node->type) {
            case NodeType::VALUE:
                str_stream << static_cast<ValueNode *>(node.get())->val;
                str += str_stream.str();
                break;
            case NodeType::VARIABLE:
                str += vars[static_cast<VariableNode *>(node.get())->idx];
                break;
            default:
                break;
        }
    } else if (!node->left != !node->right) {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        char bracketL = '(', bracketR = ')';
        switch (op) {
            case Operation::LOG:
                str += "\\log_10";
                break;
            case Operation::LN:
                str += "\\log_e";
                break;
            case Operation::ABS:
                bracketL = '|';
                bracketR = '|';
                break;
            default:
                str += "\\" + operations[static_cast<uint64_t>(op)].op_str[0];
                break;
        }
        //str += operations[static_cast<uint64_t>(op)].op_str[0];
        //if (op == Operation::MINUS) {
        //    toStringDefault(node->left ? node->left : node->right, str);
        //} else {
        str += bracketL;
        toStringLatex(node->left ? node->left : node->right, str);
        str += bracketR;
        //}
    } else {
        Operation op = static_cast<OperationNode *>(node.get())->op;
        switch (op) {
            case Operation::DIV:
                str += "\\dfrac{";
                toStringLatex(node->left, str);
                str += "}{";
                toStringLatex(node->right, str);
                str += "}";
                break;
            case Operation::POW:
                toStringLatex(node->left, str);
                str += "^{";
                toStringLatex(node->right, str);
                str += "}";
                break;
            default:
                if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
                    str += "(";
                    toStringLatex(node->left, str);
                    str += ")";
                } else {
                    toStringLatex(node->left, str);
                }
                if (op != Operation::MUL) {
                    str += " " + operations[static_cast<uint64_t>(op)].op_str[0] + " ";
                }
                if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
                    str += "(";
                    toStringLatex(node->right, str);
                    str += ")";
                } else {
                    toStringLatex(node->right, str);
                }
                break;
        }
        // if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
        //     str += "(";
        //     toStringLatex(node->left, str);
        //     str += ")";
        // } else {
        //     toStringLatex(node->left, str);
        // }
        // str += " " + operations[static_cast<uint64_t>(op)].op_str[0] + " ";
        // if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
        //     str += "(";
        //     toStringLatex(node->right, str);
        //     str += ")";
        // } else {
        //     toStringLatex(node->right, str);
        // }
    }
}

FunctionalTree::NodePtr FunctionalTree::copyTree (const NodePtr &node) const {
    if (!node) {
        return nullptr;
    }
    NodePtr copy;
    const OperationNode *operation = static_cast<const OperationNode *>(node.get());
    const ValueNode *value = static_cast<const ValueNode *>(node.get());
    const VariableNode *variable = static_cast<const VariableNode *>(node.get());
    switch (node->type) {
        case NodeType::OPERATION:
            copy = std::make_unique<OperationNode>(operation->op);
            break;
        case NodeType::VALUE:
            copy = std::make_unique<ValueNode>(value->val);
            break;
        case NodeType::VARIABLE:
            copy = std::make_unique<VariableNode>(variable->idx);
            break;
        default:
            break;
    }
    copy->priority = node->priority;
    copy->left = copyTree(node->left);
    copy->right = copyTree(node->right);
    return copy;
}

FunctionalTree::FunctionalTree (const NodePtr &node) {
    root = copyTree(node);
}

FunctionalTree::FunctionalTree (NodePtr &&tree) {
    root = std::move(tree);
}

FunctionalTree::FunctionalTree () {}

FunctionalTree::FunctionalTree (const std::string &func) {
    reset(func, {});
}

FunctionalTree::FunctionalTree (const std::string &func, const std::string &var) {
    reset(func, {var});
}

FunctionalTree::FunctionalTree (const std::string &func, const std::vector<std::string> &vars) {
    reset(func, vars);
}

FunctionalTree::FunctionalTree (const FunctionalTree &tree) {
    root = copyTree(tree.root);
    vars = tree.vars;
}

FunctionalTree::FunctionalTree (FunctionalTree &&tree) {
    root = std::move(tree.root);
    vars = std::move(tree.vars);
}

FunctionalTree::~FunctionalTree () {}

void FunctionalTree::reset (const std::string &func, const std::vector<std::string> &vars) {
    inputCheck(vars);
    this->vars = vars;
    root = buildTree(func);
}

double FunctionalTree::func (double x) const {
    return getVal(root, {x});
}

double FunctionalTree::func (const std::vector<double> &X) const {
    return getVal(root, X);
}

double FunctionalTree::calculate () const {
    return getVal(root, {});
}

std::vector<std::string> FunctionalTree::getVarList () const {
    return vars;
}

FunctionalTree FunctionalTree::getCoeff (uint64_t idx) const {
    auto tmp = root.get();
    while (tmp->left && tmp->type == NodeType::OPERATION) {
        if (tmp->right) {
            if (tmp->right->type == NodeType::VARIABLE) {
                VariableNode *var = static_cast<VariableNode *>(tmp->right.get());
                OperationNode *operation = static_cast<OperationNode *>(tmp);
                if (var->idx == idx) {
                    if (operation->op == Operation::MINUS) {
                        return FunctionalTree("-1");
                    }
                    if (operation->op == Operation::PLUS) {
                        return FunctionalTree("1");
                    }
                    if (operation->op == Operation::MUL) {
                        return FunctionalTree(tmp->left);
                    }
                }
            }
            if (tmp->right->right && tmp->right->right->type == NodeType::VARIABLE) {
                VariableNode *var = static_cast<VariableNode *>(tmp->right->right.get());
                OperationNode *operation = static_cast<OperationNode *>(tmp);
                if (var->idx == idx) {
                    FunctionalTree tr(tmp->right->left);
                    if (operation->op == Operation::MINUS) {
                        tr.root->priority = 0;
                        auto node = tr.buildTree("-");
                        tr.addToTree(tr.root, node);
                    }
                    return tr;
                }
            }
        }
        tmp = tmp->left.get();
    }
    if (tmp->type == NodeType::VARIABLE) {
        VariableNode *n = static_cast<VariableNode *>(tmp);
        if (n->idx == idx) {
            return FunctionalTree("1");
        }
    }
    return FunctionalTree("0");
}

FunctionalTree FunctionalTree::getCoeff (const std::string &param) const {
    uint64_t idx;
    auto it = std::find(vars.cbegin(), vars.cend(), param);
    idx = std::distance(vars.begin(), it);
    if (idx == vars.size()) {
        throw std::logic_error("Operation \"getCoeff\": var \"" + param + "\" not found in var list");;
    }
    return getCoeff(idx);
}

// FunctionalTree FunctionalTree::getDiv () const {
//     if (root->type == NodeType::OPERATION) {
//         return FunctionalTree(root->right);
//     } else {
//         return FunctionalTree();
//     }
// }

void FunctionalTree::printTree () const {
    printTree(root, std::cout);
}

void FunctionalTree::printFunc () const {
    printFunc(root, std::cout);
}

std::string FunctionalTree::toString (Style style) const {
    std::string str;
    if (style == Style::LATEX) {

    }
    switch (style) {
        case Style::DEFAULT:
            toStringDefault(root, str);
            break;
        case Style::GNUPLOT:
            //toStringGNUPlot(root, str);
            break;
        case Style::LATEX:
            toStringLatex(root, str);
            break;
        default:
            break;
    }
    return str;
}

//void FunctionalTree::simplify () {}

FunctionalTree &FunctionalTree::operator= (const FunctionalTree &tree) {
    if (this == &tree) {
        return *this;
    }
    root = copyTree(tree.root);
    vars = tree.vars;
    return *this;
}

FunctionalTree &FunctionalTree::operator= (FunctionalTree &&tree) {
    if (this == &tree) {
        return *this;
    }
    root = std::move(tree.root);
    vars = std::move(tree.vars);
    return *this;
}

double FunctionalTree::operator() (double x) const {
    return getVal(root, {x});
}

double FunctionalTree::operator() (const std::vector<double> &X) const {
    return getVal(root, X);
}

std::ostream &operator<< (std::ostream &output, const FunctionalTree &tree) {
    //tree.printFunc(tree.root, output);
    output << tree.toString(Style::DEFAULT);
    return output;
}

std::ifstream &operator>> (std::ifstream &file, FunctionalTree &tree) {
    std::string func;
    std::vector<std::string> vars;
    uint64_t size;
    file >> func >> size;
    for (uint64_t i = 0; i < size; ++i) {
        std::string var;
        file >> var;
        vars.push_back(var);
    }
    tree.reset(func, vars);
    return file;
}

std::ofstream &operator<< (std::ofstream &file, const FunctionalTree &tree) {
    tree.printFunc(tree.root, file);
    file << ' ' << tree.vars.size();
    for (uint64_t i = 0; i < tree.vars.size(); ++i) {
        file << ' ' << tree.vars[i];
    }
    return file;
}