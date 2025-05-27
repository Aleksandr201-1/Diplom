#include <Math/FuncMaker.hpp>

//var style

FuncMaker::OperationStruct::OperationStruct (const std::vector<std::string> &op_str, const std::function<double (double, double)> &func, uint64_t priority) 
: op_str(op_str), func(func), priority(priority) {}

FuncMaker::OperationStruct::~OperationStruct () {}

const std::vector<FuncMaker::OperationStruct> FuncMaker::operations = {
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
    {{"sign"},                     [] (double x, double y) {return x >= 0.0 ? 1 : -1;},                    1}
};

const std::map<std::string, double> FuncMaker::const_vals = {
    {"pi", std::acos(-1.0)},
    {"e", std::exp(1.0)}
};

FuncMaker::FMNode::FMNode (NodeType type) : type(type), priority(0) {}
FuncMaker::FMNode::~FMNode () {}

FuncMaker::OperationNode::OperationNode (uint64_t idx) : FMNode(NodeType::OPERATION), idx(idx) {}
FuncMaker::OperationNode::~OperationNode () {}

FuncMaker::ValueNode::ValueNode (double val) : FMNode(NodeType::VALUE), val(val) {}
FuncMaker::ValueNode::~ValueNode () {}

FuncMaker::VariableNode::VariableNode (uint64_t idx) : FMNode(NodeType::VARIABLE), idx(idx) {}
FuncMaker::VariableNode::~VariableNode () {}

void FuncMaker::inputCheck (const std::vector<std::string> &vars) const {
    for (uint64_t i = 0; i < vars.size(); ++i) {
        for (uint64_t j = 0; j < operations.size(); ++j) {
            auto it = std::find(operations[j].op_str.cbegin(), operations[j].op_str.cend(), vars[i]);
            uint64_t idx = std::distance(operations[j].op_str.cbegin(), it);
            if (idx != operations[j].op_str.size()) {
                throw std::runtime_error("Operation \"inputCheck\": the variable name \"" + vars[i] + "\" is used, which is a name of operation");
            }
        }
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        auto it = const_vals.find(vars[i]);
        if (it != const_vals.end()) {
            throw std::runtime_error("Operation \"inputCheck\": the variable name \"" + it->first + "\" is used, which is a name of constant");
        }
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        auto it = vals.find(vars[i]);
        if (it != vals.end()) {
            throw std::runtime_error("Operation \"inputCheck\": the variable name \"" + it->first + "\" is used, which is a name of value");
        }
    }
    for (uint64_t i = 0; i < vars.size(); ++i) {
        uint64_t count = std::count(vars.cbegin() + i, vars.cend(), vars[i]);
        if (count != 1) {
            throw std::runtime_error("Operation \"inputCheck\": the variable name \"" + vars[i] + "\" is used multiple times");
        }
    }
}

std::string FuncMaker::readOperation (const std::string &func, uint64_t &i) const {
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

std::string FuncMaker::readWord (const std::string &func, uint64_t &i) const {
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

double FuncMaker::readNumber (const std::string &func, uint64_t &i) const {
    std::string str;
    while (((func[i] >= '0' && func[i] <= '9') || func[i] == '.') && i < func.size()) {
        str += func[i];
        ++i;
    }
    return std::atof(str.c_str());
}

std::string FuncMaker::readInbrace (const std::string &func, uint64_t &i) const {
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

uint64_t FuncMaker::getOperation (const std::string &str) const {
    uint64_t op = -1;//Operation::NOT_AN_OPERATION;
    //std::cout << "op size: " << operations.size() << "\n";
    //std::cout << "[lus]: " << operations[0].op_str[0] << "\n";
    for (uint64_t i = 0; i < operations.size(); ++i) {
        auto it = std::find(operations[i].op_str.cbegin(), operations[i].op_str.cend(), str);
        uint64_t idx = std::distance(operations[i].op_str.cbegin(), it);
        //std::cout << "idx = " << idx << "\n";
        //std::cout << str << "==" << operations[i].op_str[idx] << "\n";
        if (idx != operations[i].op_str.size()) {
            op = i;
            //return idx;
            break;
            //throw std::runtime_error("Operation \"inputCheck\": the variable name \"" + vars[i] + "\" is used, which is a name of operation");
        }
        // for (uint64_t j = 0; j < operations[i].op_str.size(); ++j) {
        //     if (str == operations[i].op_str[j]) {
        //         op = static_cast<Operation>(i);
        //         break;
        //     }
        // }
    }
    return op;
}

double FuncMaker::getConstant (const std::string &str) const {
    double constant = std::nan("1");
    auto it = const_vals.find(str);
    if (it != const_vals.end()) {
        constant = it->second;
    }
    return constant;
}

uint64_t FuncMaker::getVariable (const std::string &str) const {
    auto it = std::find(vars.cbegin(), vars.cend(), str);
    return std::distance(vars.cbegin(), it);
}

uint64_t FuncMaker::getPriority (uint64_t idx) const {
    //uint64_t idx = static_cast<uint64_t>(op);
    return operations[idx].priority;
}

double FuncMaker::useOperation (uint64_t idx, double x, double y) const {
    //uint64_t idx = static_cast<uint64_t>(op);
    return operations[idx].func(x, y);
}

double FuncMaker::calcNode (const NodePtr &node, const std::vector<double> &X) const {
    if (!node->left.get() && !node->right.get()) {
        if (node->type == NodeType::VALUE) {
            return static_cast<ValueNode *>(node.get())->val;
        } else {
            return X[static_cast<VariableNode *>(node.get())->idx];
        }
    }
    OperationNode *operation = static_cast<OperationNode *>(node.get());
    if (node->left && !node->right) {
        return useOperation(operation->idx, calcNode(operation->left, X), 0);
    }
    if (!node->left && node->right) {
        return useOperation(operation->idx, calcNode(operation->right, X), 0);
    }
    double a = calcNode(node->left, X);
    double b = calcNode(node->right, X);
    return useOperation(operation->idx, a, b);
}

void FuncMaker::addToTree (NodePtr &tree, NodePtr &node) {
    if (!tree.get()) {
        if (node->type == NodeType::OPERATION && node->priority != 0) {
            OperationNode *operation = static_cast<OperationNode *>(node.get());
            if (operations[operation->idx].op_str[0] == "-") {
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

FuncMaker::NodePtr FuncMaker::buildTree (const std::string &func) {
    std::string tmp;
    uint64_t i = 0;

    double num;
    uint64_t op_idx;
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
            op_idx = getOperation(tmp);
            //std::cout << "op: " << op_idx << "\n";
            if (op_idx == uint64_t(-1)) {
                double constant = getConstant(tmp);
                if (std::isnan(constant)) {
                    constant = getValue(tmp);
                }
                if (std::isnan(constant)) {
                    idx = getVariable(tmp);
                    if (idx == vars.size()) {
                        throw std::runtime_error("Operation \"buildTree\": variable \"" + tmp + "\" not found in var list");
                    }
                    current = std::make_unique<VariableNode>(idx);
                } else {
                    current = std::make_unique<ValueNode>(constant);
                }
            } else {
                current = std::make_unique<OperationNode>(op_idx);
                current->priority = getPriority(op_idx);
            }
        }
        addToTree(node, current);
    }
    return node;
}

void FuncMaker::printTree (const NodePtr &node, std::ostream &out) const {
    if (!node) {
        return;
    }
    //uint64_t op_idx;
    switch (node->type) {
        case NodeType::VALUE:
            out << static_cast<ValueNode *>(node.get())->val;
            break;
        case NodeType::VARIABLE:
            out << vars[static_cast<VariableNode *>(node.get())->idx];
            break;
        case NodeType::OPERATION:
            //op_idx = static_cast<OperationNode *>(node.get())->idx;
            out << static_cast<OperationNode *>(node.get())->idx;
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

void FuncMaker::printFunc (const NodePtr &node, std::ostream &out) const {
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
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        out << operations[op_idx].op_str[0];
        if (operations[op_idx].op_str[0] == "-") {
            printFunc(node->left ? node->left : node->right, out);
        } else {
            out << "(";
            printFunc(node->left ? node->left : node->right, out);
            out << ")";
        }
    } else {
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
            out << "(";
            printFunc(node->left, out);
            out << ")";
        } else {
            printFunc(node->left, out);
        }
        out << operations[op_idx].op_str[0];
        if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
            out << "(";
            printFunc(node->right, out);
            out << ")";
        } else {
            printFunc(node->right, out);
        }
    }
}

void FuncMaker::toStringDefault (const NodePtr &node, std::string &str) const {
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
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        str += operations[op_idx].op_str[0];
        if (operations[op_idx].op_str[0] == "-") {
            toStringDefault(node->left ? node->left : node->right, str);
        } else {
            str += "(";
            toStringDefault(node->left ? node->left : node->right, str);
            str += ")";
        }
    } else {
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
            str += "(";
            toStringDefault(node->left, str);
            str += ")";
        } else {
            toStringDefault(node->left, str);
        }
        str += " " + operations[op_idx].op_str[0] + " ";
        if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
            str += "(";
            toStringDefault(node->right, str);
            str += ")";
        } else {
            toStringDefault(node->right, str);
        }
    }
}

//void FuncMaker::toStringGNUPlot (const NodePtr &node, std::string &str) const {}

void FuncMaker::toStringLatex (const NodePtr &node, std::string &str) const {
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
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        std::string op_str = operations[op_idx].op_str[0];
        char bracketL = '(', bracketR = ')';
        if (op_str == "log") {
            str += "\\log_10";
        } else if (op_str == "ln") {
            str += "\\log_e";
        } else if (op_str == "abs") {
            bracketL = '|';
            bracketR = '|';
        } else {
            str += "\\" + operations[op_idx].op_str[0];
        }
        // switch (op_idx) {
        //     case Operation::LOG:
        //         str += "\\log_10";
        //         break;
        //     case Operation::LN:
        //         str += "\\log_e";
        //         break;
        //     case Operation::ABS:
        //         bracketL = '|';
        //         bracketR = '|';
        //         break;
        //     default:
        //         str += "\\" + operations[op_idx].op_str[0];
        //         break;
        // }
        //======================================
        //str += operations[op_idx].op_str[0];
        //if (op == Operation::MINUS) {
        //    toStringDefault(node->left ? node->left : node->right, str);
        //} else {
        str += bracketL;
        toStringLatex(node->left ? node->left : node->right, str);
        str += bracketR;
        //}
    } else {
        uint64_t op_idx = static_cast<OperationNode *>(node.get())->idx;
        std::string op_str = operations[op_idx].op_str[0];
        if (op_str == "/") {
            str += "\\dfrac{";
            toStringLatex(node->left, str);
            str += "}{";
            toStringLatex(node->right, str);
            str += "}";
        } else if (op_str == "^") {
            toStringLatex(node->left, str);
            str += "^{";
            toStringLatex(node->right, str);
            str += "}";
        } else {
            if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
                str += "(";
                toStringLatex(node->left, str);
                str += ")";
            } else {
                toStringLatex(node->left, str);
            }
            if (op_str != "*") {
                str += " " + operations[op_idx].op_str[0] + " ";
            }
            if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
                str += "(";
                toStringLatex(node->right, str);
                str += ")";
            } else {
                toStringLatex(node->right, str);
            }
        }
        // switch (op) {
        //     case Operation::DIV:
        //         str += "\\dfrac{";
        //         toStringLatex(node->left, str);
        //         str += "}{";
        //         toStringLatex(node->right, str);
        //         str += "}";
        //         break;
        //     case Operation::POW:
        //         toStringLatex(node->left, str);
        //         str += "^{";
        //         toStringLatex(node->right, str);
        //         str += "}";
        //         break;
        //     default:
        //         if (node->left->priority == 0 && node->left->type == NodeType::OPERATION) {
        //             str += "(";
        //             toStringLatex(node->left, str);
        //             str += ")";
        //         } else {
        //             toStringLatex(node->left, str);
        //         }
        //         if (op != Operation::MUL) {
        //             str += " " + operations[static_cast<uint64_t>(op)].op_str[0] + " ";
        //         }
        //         if (node->right->priority == 0 && node->right->type == NodeType::OPERATION) {
        //             str += "(";
        //             toStringLatex(node->right, str);
        //             str += ")";
        //         } else {
        //             toStringLatex(node->right, str);
        //         }
        //         break;
        // }
    }
}

FuncMaker::NodePtr FuncMaker::copyTree (const NodePtr &node) const {
    if (!node) {
        return nullptr;
    }
    NodePtr copy;
    const OperationNode *operation = static_cast<const OperationNode *>(node.get());
    const ValueNode *value = static_cast<const ValueNode *>(node.get());
    const VariableNode *variable = static_cast<const VariableNode *>(node.get());
    switch (node->type) {
        case NodeType::OPERATION:
            copy = std::make_unique<OperationNode>(operation->idx);
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

FuncMaker::FuncMaker (const NodePtr &node) {
    root = copyTree(node);
}

FuncMaker::FuncMaker (NodePtr &&tree) {
    root = std::move(tree);
}

FuncMaker::FuncMaker () {}

FuncMaker::FuncMaker (const std::string &func) {
    reset(func, {});
}

FuncMaker::FuncMaker (const std::string &func, const std::string &var) {
    reset(func, {var});
}

FuncMaker::FuncMaker (const std::string &func, const std::vector<std::string> &vars) {
    reset(func, vars);
}

FuncMaker::FuncMaker (const std::string &func, const std::initializer_list<std::string> &vars) {
    reset(func, vars);
}

FuncMaker::FuncMaker (const FuncMaker &tree) {
    root = copyTree(tree.root);
    vars = tree.vars;
}

FuncMaker::FuncMaker (FuncMaker &&tree) {
    root = std::move(tree.root);
    vars = std::move(tree.vars);
}

FuncMaker::~FuncMaker () {}

void FuncMaker::reset (const std::string &func, const std::initializer_list<std::string> &vars) {
    reset(func, std::vector<std::string>(vars));
}

void FuncMaker::reset (const std::string &func, const std::vector<std::string> &vars) {
    inputCheck(vars);
    this->vars = vars;
    root = buildTree(func);
}

void FuncMaker::setValue (const std::string &name, double val) {
    for (uint64_t j = 0; j < operations.size(); ++j) {
        auto it = std::find(operations[j].op_str.cbegin(), operations[j].op_str.cend(), name);
        uint64_t idx = std::distance(operations[j].op_str.cbegin(), it);
        if (idx != operations[j].op_str.size()) {
            throw std::runtime_error("Operation \"setValue\": the value name \"" + name + "\" is used, which is a name of operation");
        }
    }
    if (const_vals.find(name) != const_vals.end()) {
        throw std::runtime_error("Operation \"setValue\": the value name \"" + name + "\" is used, which is a name of constant");
    }
    if (vals.find(name) != vals.end()) {
        throw std::runtime_error("Operation \"setValue\": the value name \"" + name + "\" is used multiple times");
    }
    auto it = std::find(vars.cbegin(), vars.cend(), name);
    if (it != vars.end()) {
        throw std::runtime_error("Operation \"setValue\": the value name \"" + name + "\" is used, which is a name of variable");
    }
    vals[name] = val;
}

double FuncMaker::getValue (const std::string &name) const {
    double ans = std::nan("1");
    auto it = vals.find(name);
    if (it != vals.end()) {
        ans = it->second;
    }
    return ans;
}

double FuncMaker::func (double x) const {
    return calcNode(root, {x});
}

double FuncMaker::func (const std::vector<double> &X) const {
    return calcNode(root, X);
}

double FuncMaker::func (const std::initializer_list<double> &X) const {
    return calcNode(root, X);
}

double FuncMaker::calculate () const {
    return calcNode(root, {});
}

std::vector<std::string> FuncMaker::getVariableList () const {
    return vars;
}

std::vector<std::string> FuncMaker::getValueList () const {
    std::vector<std::string> vals_str;
    for (const auto &it : vals) {
        vals_str.push_back(it.first);
    }
    return vals_str;
}

FuncMaker FuncMaker::getCoeff (uint64_t idx) const {
    auto tmp = root.get();
    while (tmp->left && tmp->type == NodeType::OPERATION) {
        if (tmp->right) {
            if (tmp->right->type == NodeType::VARIABLE) {
                VariableNode *var = static_cast<VariableNode *>(tmp->right.get());
                OperationNode *operation = static_cast<OperationNode *>(tmp);
                uint64_t op_idx = operation->idx;
                std::string op_str = operations[op_idx].op_str[0];
                if (var->idx == idx) {
                    if (op_str == "-") {
                        return FuncMaker("-1");
                    }
                    if (op_str == "+") {
                        return FuncMaker("1");
                    }
                    if (op_str == "*") {
                        return FuncMaker(tmp->left);
                    }
                }
            }
            if (tmp->right->right && tmp->right->right->type == NodeType::VARIABLE) {
                VariableNode *var = static_cast<VariableNode *>(tmp->right->right.get());
                OperationNode *operation = static_cast<OperationNode *>(tmp);
                if (var->idx == idx) {
                    FuncMaker tr(tmp->right->left);
                    uint64_t op_idx = operation->idx;
                    std::string op_str = operations[op_idx].op_str[0];
                    if (op_str == "-") {
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
            return FuncMaker("1");
        }
    }
    return FuncMaker("0");
}

FuncMaker FuncMaker::getCoeff (const std::string &param) const {
    uint64_t idx;
    auto it = std::find(vars.cbegin(), vars.cend(), param);
    idx = std::distance(vars.begin(), it);
    if (idx == vars.size()) {
        throw std::runtime_error("Operation \"getCoeff\": var \"" + param + "\" not found in var list");
    }
    return getCoeff(idx);
}

// FuncMaker FuncMaker::getDiv () const {
//     if (root->type == NodeType::OPERATION) {
//         return FuncMaker(root->right);
//     } else {
//         return FuncMaker();
//     }
// }

void FuncMaker::printTree () const {
    printTree(root, std::cout);
}

void FuncMaker::printFunc () const {
    printFunc(root, std::cout);
}

std::string FuncMaker::toString (Style style) const {
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

//void FuncMaker::simplify () {}

FuncMaker &FuncMaker::operator= (const FuncMaker &tree) {
    if (this == &tree) {
        return *this;
    }
    root = copyTree(tree.root);
    vars = tree.vars;
    return *this;
}

FuncMaker &FuncMaker::operator= (FuncMaker &&tree) {
    if (this == &tree) {
        return *this;
    }
    root = std::move(tree.root);
    vars = std::move(tree.vars);
    return *this;
}

double FuncMaker::operator() (double x) const {
    return calcNode(root, {x});
}

double FuncMaker::operator() (const std::vector<double> &X) const {
    return calcNode(root, X);
}

double FuncMaker::operator() (const std::initializer_list<double> &X) const {
    return calcNode(root, X);
}

std::istream &operator>> (std::istream &input, FuncMaker &tree) {
    std::string str;
    input >> str;
    return input;
}

std::ostream &operator<< (std::ostream &output, const FuncMaker &tree) {
    //tree.printFunc(tree.root, output);
    output << tree.toString(FuncMaker::Style::DEFAULT);
    return output;
}

std::ifstream &operator>> (std::ifstream &file, FuncMaker &tree) {
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

std::ofstream &operator<< (std::ofstream &file, const FuncMaker &tree) {
    tree.printFunc(tree.root, file);
    file << ' ' << tree.vars.size();
    for (uint64_t i = 0; i < tree.vars.size(); ++i) {
        file << ' ' << tree.vars[i];
    }
    return file;
}