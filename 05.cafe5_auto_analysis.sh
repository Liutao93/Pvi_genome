#!/bin/bash

# CAFE5 自动分析脚本 (从无数字后缀的初始文件开始)
# 功能：自动识别并复制初始文件，交互式配置参数，循环运行直至分析完成

set -euo pipefail  # 更严格的错误处理

# 初始化变量
CONFIG_FILE="cafe5_analysis.config"
SCRIPT_VERSION="2.0"
ORIGINAL_FILE="Orthogroups.GeneCount.cafe.filtered.tsv"
STARTING_FILE="Orthogroups.GeneCount.cafe.filtered0.tsv"

echo "=================================================="
echo "    CAFE5 自动分析脚本 v$SCRIPT_VERSION"
echo "=================================================="
echo "目标：从 $ORIGINAL_FILE 开始，自动迭代直至分析完成"
echo "=================================================="

# 函数：检查初始文件并准备
prepare_initial_file() {
    echo "检查初始文件..."
    
    if [ ! -f "$ORIGINAL_FILE" ]; then
        echo " 错误：未找到初始文件 $ORIGINAL_FILE"
        echo " 请确保该文件存在于当前目录"
        return 1
    fi
    
    # 如果起始文件不存在，从原始文件创建
    if [ ! -f "$STARTING_FILE" ]; then
        echo "创建起始文件: $ORIGINAL_FILE → $STARTING_FILE"
        cp "$ORIGINAL_FILE" "$STARTING_FILE"
        
        # 验证复制是否成功
        if [ $? -eq 0 ] && [ -s "$STARTING_FILE" ]; then
            echo " 起始文件创建成功"
            return 0
        else
            echo " 起始文件创建失败"
            return 1
        fi
    else
        echo "  起始文件已存在: $STARTING_FILE"
        return 0
    fi
}

# 函数：交互式配置
setup_parameters() {
    echo ""
    echo "=== 参数配置 ==="
    
    # 树文件配置
    read -p "请输入树文件路径 [默认: cafe_formated.nwk]: " TREE_FILE
    TREE_FILE=${TREE_FILE:-cafe_formated.nwk}
    
    if [ ! -f "$TREE_FILE" ]; then
        echo " 错误：树文件 $TREE_FILE 不存在"
        return 1
    fi
    
    # 输出前缀
    read -p "请输入输出前缀 [默认: cafe_result]: " OUTPUT_PREFIX
    OUTPUT_PREFIX=${OUTPUT_PREFIX:-cafe_result}
    
    # 日志文件
    read -p "请输入日志文件名 [默认: cafe_analysis.log]: " LOG_FILE
    LOG_FILE=${LOG_FILE:-cafe_analysis.log}
    
    # 核心数
    read -p "请输入使用的CPU核心数 [默认: 20]: " CORES
    CORES=${CORES:-20}
    
    # p值
    read -p "请输入p值阈值 [默认: 0.05]: " PVALUE
    PVALUE=${PVALUE:-0.05}
    
    # 初始文件编号
    CURRENT_FILTER_NUM=0
    ITERATION_COUNT=0
    
    echo ""
    echo "=== 配置总结 ==="
    echo " 起始文件: $STARTING_FILE"
    echo " 树文件: $TREE_FILE"
    echo " 输出前缀: $OUTPUT_PREFIX"
    echo " 日志文件: $LOG_FILE"
    echo "  CPU核心: $CORES"
    echo " p值阈值: $PVALUE"
    echo ""
    
    read -p "确认开始分析? [y/N]: " confirm
    case "$confirm" in
        [yY]|[yY][eE][sS])
            echo " 开始分析..."
            return 0
            ;;
        *)
            echo "分析已取消"
            return 1
            ;;
    esac
}

# 函数：检查CAFE5是否可用
check_cafe5() {
    if ! command -v cafe5 &> /dev/null; then
        echo " 错误：未找到 cafe5 命令"
        echo " 请确保CAFE5已正确安装并在PATH中"
        return 1
    fi
    echo " CAFE5 可用"
    return 0
}

# 函数：运行单次迭代
run_iteration() {
    local current_num=$1
    local iteration=$2
    
    echo ""
    echo "第 $iteration 次迭代"
    
    # 检查是否已完成
    if [ -d "$OUTPUT_PREFIX" ] && [ -n "$(ls -A "$OUTPUT_PREFIX" 2>/dev/null)" ]; then
        echo "分析完成！$OUTPUT_PREFIX 文件夹不为空"
        return 2  # 特殊返回码表示完成
    fi
    
    # 设置输入输出文件
    local input_file="Orthogroups.GeneCount.cafe.filtered${current_num}.tsv"
    local next_num=$((current_num + 1))
    local output_file="Orthogroups.GeneCount.cafe.filtered${next_num}.tsv"
    
    echo " 输入文件: $input_file"
    echo " 输出文件: $output_file"
    
    # 检查输入文件是否存在
    if [ ! -f "$input_file" ]; then
        echo " 错误：输入文件 $input_file 不存在"
        return 1
    fi
    
    # 如果有日志文件，提取OG行
    if [ -f "$LOG_FILE" ]; then
        echo "从日志文件提取OG家族..."
        
        # 提取OG行并显示
        local og_count=$(grep -c "^OG" "$LOG_FILE" 2>/dev/null || true)
        
        if [ $og_count -gt 0 ]; then
            echo "发现 $og_count 个需要过滤的OG家族:"
            grep "^OG" "$LOG_FILE" | awk -F":" '{print $1}' | tee "family.remove.id.txt"
            echo " 家族ID已保存到: family.remove.id.txt"
        else
            echo "  日志文件中未找到OG家族行"
            # 创建空文件
            > "family.remove.id.txt"
        fi
    else
        echo "  日志文件不存在（可能是首次运行）"
        > "family.remove.id.txt"
    fi
    
    # 执行过滤（只有当family.remove.id.txt非空时才过滤）
    if [ -s "family.remove.id.txt" ]; then
        echo " 执行基因家族过滤..."
        grep -v -f "family.remove.id.txt" "$input_file" > "$output_file"
        
        # 检查过滤结果
        if [ ! -s "$output_file" ]; then
            echo " 错误：过滤后的文件为空"
            return 1
        fi
        
        local input_lines=$(wc -l < "$input_file")
        local output_lines=$(wc -l < "$output_file")
        local removed_lines=$((input_lines - output_lines))
        
        echo " 过滤完成: 移除 $removed_lines 行，剩余 $output_lines 行"
    else
        echo "  无家族需要过滤，直接复制文件"
        cp "$input_file" "$output_file"
    fi
    
    # 运行CAFE5分析
    echo " 运行CAFE5分析..."
    local cafe_cmd="cafe5 --infile \"$output_file\" --tree \"$TREE_FILE\" --output_prefix \"$OUTPUT_PREFIX\" --cores \"$CORES\" -p --pvalue \"$PVALUE\""
    echo "执行命令: $cafe_cmd &> \"$LOG_FILE\""
    
    eval $cafe_cmd &> "$LOG_FILE"
    local cafe_exit=$?
    
    if [ $cafe_exit -ne 0 ]; then
        echo " CAFE5退出代码: $cafe_exit (可能正常)"
    else
        echo " CAFE5分析完成"
    fi
    
    # 返回下一个文件编号
    echo $next_num
    return 0
}

# 主函数
main() {
    echo "开始时间: $(date)"
    echo "当前目录: $(pwd)"
    echo ""
    
    # 检查初始文件
    if ! prepare_initial_file; then
        exit 1
    fi
    
    # 检查CAFE5
    if ! check_cafe5; then
        exit 1
    fi
    
    # 交互式配置
    if ! setup_parameters; then
        exit 0
    fi
    
    # 初始化变量
    local current_filter_num=0
    local iteration_count=0
    local max_safety_iterations=1000  # 安全上限
    
    echo ""
    echo "=== 开始主循环 ==="
    
    # 主循环
    while [ $iteration_count -lt $max_safety_iterations ]; do
        iteration_count=$((iteration_count + 1))
        
        # 运行单次迭代
        if ! result=$(run_iteration $current_filter_num $iteration_count); then
            case $? in
                1)
                    echo " 迭代失败，停止运行"
                    break
                    ;;
                2)
                    echo " 分析成功完成！"
                    break
                    ;;
            esac
        else
            current_filter_num=$result
            echo " 迭代 $iteration_count 完成，当前文件编号: $current_filter_num"
        fi
        
        echo "----------------------------------------"
    done
    
    # 最终总结
    echo ""
    echo "=================================================="
    echo "                分析完成总结"
    echo "=================================================="
    echo " 总迭代次数: $iteration_count"
    echo " 最终过滤文件: Orthogroups.GeneCount.cafe.filtered${current_filter_num}.tsv"
    echo " 输出目录: $OUTPUT_PREFIX"
    echo " 日志文件: $LOG_FILE"
    echo "完成时间: $(date)"
    
    if [ $iteration_count -eq $max_safety_iterations ]; then
        echo "  达到安全迭代上限，请检查分析状态"
    fi
    
    # 显示最终结果
    if [ -d "$OUTPUT_PREFIX" ] && [ -n "$(ls -A "$OUTPUT_PREFIX")" ]; then
        echo ""
        echo " 输出目录内容:"
        ls -la "$OUTPUT_PREFIX"
    fi
}

# 信号处理
trap 'echo ""; echo "脚本被用户中断"; exit 1' INT TERM

# 运行主程序
main
