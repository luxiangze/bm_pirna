import os
import json
import glob
import sys
import argparse

def generate_group_files(bam_dir, config_file, output_dir):
    # 读取配置文件
    print(f"正在读取配置文件: {config_file}")
    with open(config_file, 'r') as f:
        config = json.load(f)
    sample_groups = config.get('samples', {})
    print(f"配置中的样本组: {sample_groups}")
    
    # 查找BAM文件 - 使用更精确的匹配方式
    bam_files = []
    for sample_dir in os.listdir(bam_dir):
        bam_path = os.path.join(bam_dir, sample_dir, "Aligned.sortedByCoord.out.bam")
        if os.path.isfile(bam_path):
            bam_files.append(bam_path)
            print(f"找到BAM文件: {bam_path}")
    
    print(f"总共找到 {len(bam_files)} 个BAM文件")
    
    os.makedirs(output_dir, exist_ok=True)
    for group_name, samples in sample_groups.items():
        group_bams = []
        for bam in bam_files:
            # 获取样本目录名称
            sample_dir = os.path.basename(os.path.dirname(bam))
            print(f"检查BAM文件 {bam} (样本目录: {sample_dir})")
            for sample in samples:
                if sample in sample_dir:
                    group_bams.append(bam)
                    print(f"将样本 {sample_dir} 分配到组 {group_name}")
                    break
        
        if not group_bams:
            print(f"警告: 组 {group_name} 没有找到匹配的BAM文件")
            # 创建一个空文件，防止工作流失败
            group_file = os.path.join(output_dir, f"{group_name}.txt")
            with open(group_file, 'w') as f:
                f.write('')
            print(f"为组 {group_name} 创建了空文件: {group_file}")
        else:
            group_file = os.path.join(output_dir, f"{group_name}.txt")
            with open(group_file, 'w') as f:
                f.write(','.join(group_bams))
            print(f"为组 {group_name} 创建了文件: {group_file} 包含 {len(group_bams)} 个样本")
            print(f"文件内容: {','.join(group_bams)}")
        
        # 确认文件已创建
        if os.path.exists(group_file):
            print(f"确认文件存在: {group_file}")
        else:
            print(f"错误: 文件未创建: {group_file}")
            
    return True

def generate_comparison_files(config_file, output_dir):
    print(f"正在读取比较组配置: {config_file}")
    with open(config_file, 'r') as f:
        config = json.load(f)
    comparisons = config.get('comparisons', [])
    print(f"配置中定义的比较组: {comparisons}")
    
    comparison_list = []
    os.makedirs(output_dir, exist_ok=True)
    
    # 创建每个比较组的文件
    for i, (group1, group2) in enumerate(comparisons):
        comparison_name = f"comparison_{i+1}_{group1}_vs_{group2}"
        comparison_file = os.path.join(output_dir, f"{comparison_name}.txt")
        with open(comparison_file, 'w') as f:
            f.write(f"{group1}\t{group2}")
        
        # 确认文件创建成功
        if os.path.exists(comparison_file):
            print(f"成功创建比较组文件: {comparison_file}")
        else:
            print(f"错误: 未能创建比较组文件: {comparison_file}")
            
        comparison_list.append(comparison_name)
        print(f"添加比较组: {comparison_name}: {group1} vs {group2}")
    
    # 创建比较索引文件
    comparison_index = os.path.join(output_dir, "comparison_index.txt")
    try:
        with open(comparison_index, 'w') as f:
            for i, (group1, group2) in enumerate(comparisons):
                comparison_name = f"comparison_{i+1}_{group1}_vs_{group2}"
                f.write(f"{comparison_name}\t{group1}\t{group2}\n")
        
        # 确认索引文件创建成功
        if os.path.exists(comparison_index):
            print(f"成功创建比较索引文件: {comparison_index}")
            with open(comparison_index, 'r') as f:
                content = f.read()
            print(f"索引文件内容: \n{content}")
        else:
            print(f"错误: 未能创建比较索引文件: {comparison_index}")
    except Exception as e:
        print(f"创建比较索引文件时出错: {e}")
    
    return comparison_list

def parse_args():
    parser = argparse.ArgumentParser(description='生成样本组文件和比较组索引')
    parser.add_argument('--bam-dir', default='star', help='BAM文件目录')
    parser.add_argument('--config', default='config.json', help='配置文件路径')
    parser.add_argument('--output-dir', default='groups', help='输出目录')
    return parser.parse_args()

def main():
    try:
        # 将所有输出写入标准输出和标准错误输出
        sys.stdout.write("====== 开始生成组文件 ======\n")
        sys.stdout.flush()
        
        # 解析命令行参数
        args = parse_args()
        sys.stdout.write(f"输入参数: BAM目录={args.bam_dir}, 配置文件={args.config}, 输出目录={args.output_dir}\n")
        
        # 确保输出目录存在
        os.makedirs(args.output_dir, exist_ok=True)
        sys.stdout.write(f"已创建输出目录: {args.output_dir}\n")
        
        # 生成样本组文件
        sys.stdout.write("\n---- 生成样本组文件 ----\n")
        generate_group_files(args.bam_dir, args.config, args.output_dir)
        
        # 生成比较组索引文件
        sys.stdout.write("\n---- 生成比较组索引文件 ----\n")
        generate_comparison_files(args.config, args.output_dir)
        
        # 检查生成的文件
        sys.stdout.write("\n---- 检查生成的文件 ----\n")
        if os.path.exists(args.output_dir):
            files = os.listdir(args.output_dir)
            sys.stdout.write(f"在 {args.output_dir} 目录下找到的文件: {files}\n")
        else:
            sys.stdout.write(f"错误: 目录 {args.output_dir} 不存在!\n")
            
        sys.stdout.write("====== 完成生成组文件 ======\n")
        sys.stdout.flush()
    except Exception as e:
        sys.stderr.write(f"\n!! 脚本执行遇到错误: {e} !!\n")
        import traceback
        sys.stderr.write(traceback.format_exc())
        sys.stderr.flush()

if __name__ == "__main__":
    main()