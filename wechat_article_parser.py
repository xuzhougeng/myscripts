#!/usr/bin/env python

"""
provide the following parameters

- Cookie
- user_agent
- fakeid
- token

"""

import json
import requests
import time
import random
import os
import yaml
import sys
import argparse
from bs4 import BeautifulSoup

def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-j", "--json", help="json file", required=False)
    group.add_argument("-y", "--yaml", help="configuration file", required=True)

    return group.parse_args()

# 毫秒数转日期
def get_date(times):
    # print(times)
    timearr = time.localtime(times)
    date = time.strftime("%Y-%m-%d %H:%M:%S", timearr)
    return date

# 保存结果为JSON
def write_url(app_msg_list, json_name):
    with open(json_name, "w") as file:
        file.write(json.dumps(app_msg_list, indent=2, ensure_ascii=False))

# get the configure from YAML
"""
  - cookie
  - user_agent
From /mp/getappmsgext
  - key
  - pass_ticket
  - appmsg_token
"""
def read_config(yaml_file):
    with open(yaml_file, "r") as file:
        file_data = file.read()
    config = yaml.safe_load(file_data) 
    return config

# get the artile url 
def get_article_url(config, json_file):
    # set header
    headers = {
        "Cookie": config['cookie'],
        "User-Agent": config['user_agent'] 
    }

    # set the request paramters
    url = "https://mp.weixin.qq.com/cgi-bin/appmsg"
    begin = "0"
    params = {
        "action": "list_ex",
        "begin": begin,
        "count": "5",
        "fakeid": config['fakeid'],
        "type": "9",
        "token": config['token'],
        "lang": "zh_CN",
        "f": "json",
        "ajax": "1"
    }

    # read the exists result
    if os.path.exists(json_file):
        with open(json_file, "r") as file:
            app_msg_list = json.load(file)
    else:
        app_msg_list = []

    # 在不知道公众号有多少文章的情况下，使用while语句
    # 也方便重新运行时设置页数
    i = len(app_msg_list) 
    while True:
        begin = i * 5
        params["begin"] = str(begin)
        print(i)
        # 随机暂停几秒，避免过快的请求导致过快的被查到
        time.sleep(random.randint(1,10))
        resp = requests.get(url, headers=headers, params = params, verify=False)
        # 微信流量控制, 退出
        if resp.json()['base_resp']['ret'] == 200013:
            print("frequencey control, stop at {}".format(str(begin)))
            break
        #print(resp.json())
        
        # 如果返回的内容中为空则结束
        if len(resp.json()['app_msg_list']) == 0:
            print("all ariticle parsed")
            break
            
        app_msg_list.append(resp.json())
        # 翻页
        i += 1
    return app_msg_list

# get the aritle 
def get_one_article(url, config):
    # set header
    headers = {
        "Cookie": config['cookie'],
        "User-Agent": config['user_agent'] 
    }
    resp = requests.get(url, headers = headers)
    soup = BeautifulSoup(resp.text, 'lxml')

    # 获取 name='author'的tag
    author_tag = soup.find(attrs={'name':'author'}) 
    # 提取content的内容
    author = author_tag.attrs['content']
    # artile
    article_tag = soup.find(attrs={'property':'og:title'}) 
    article = article_tag.attrs['content']

    content = f'{article}\t{author}'
    print(content)
    return content
	
def get_all_article(json_file, config):

    #根据之前获取的Link进行龟速爬取
    app_msg_list = json.load(open(json_file))
    # 保存的数据列表
    info_list = []
    # 由于我们静态保存了所有链接，因此重新运行时只需要跳过前N个
    resume_site = len(info_list)
    count = 0
    for msg in app_msg_list:
        if "app_msg_list" in msg:
            for item in msg["app_msg_list"]:
                if count < resume_site:
                    count += 1
                    continue
                time.sleep(random.randint(1,4))
                link = item['link']
                #aid  = item['aid']    #文章标识符
                #title = item['title'] 
                content = get_one_article(link, config)
                info_list.append(content)
    return info_list

# get the read number, likes and so on
def get_one_appmsgstat(link, headers, key, pass_ticket, appmsg_token):

    url = "http://mp.weixin.qq.com/mp/getappmsgext"
    # POST信息
    data = {
        "is_only_read": "1",
        "is_temp_url": "0",
        "reward_uin_count": "0"
    }
    
    # 参数信息
    mid = link.split("&")[1].split("=")[1]
    idx = link.split("&")[2].split("=")[1]
    sn = link.split("&")[3].split("=")[1]
    _biz = link.split("&")[0].split("_biz=")[1]
    uin = "MjcxNDIwMDk0MA=="
    params = {
        "__biz": _biz,
        "mid": mid,
        "sn": sn,
        "idx": idx,
        "key": key, # 20-30分钟变动一次
        "pass_ticket": pass_ticket,
        "appmsg_token": appmsg_token,
        "uin": uin,
        "wxtoken": "777",
        }
    #print(params)
    content = requests.post(url, headers=headers, data=data, params=params).json()
    return content

def get_all_appmsgstat(json_file, config):
 
    # set header
    headers = {
        "Cookie": config['cookie'],
        "User-Agent": config['user_agent']
    }
    # 注意key很容易会过期, 需要及时更新
    key = config['key']
    pass_ticket = config['pass_ticket']
    appmsg_token = config['appmsg_token']

    #根据之前获取的Link进行龟速爬取
    app_msg_list = json.load(open(json_file))
    # 保存的数据列表
    info_list = []
    # 由于我们静态保存了所有链接，因此重新运行时只需要跳过前N个
    resume_site = len(info_list)
    count = 0
    for msg in app_msg_list:
        if "app_msg_list" in msg:
            for item in msg["app_msg_list"]:
                if count < resume_site:
                    count += 1
                    continue
                time.sleep(random.randint(3,7))
                aid  = item['aid']    #文章标识符
                title = item['title'] 
                print("Process {}:{}".format(aid, title))
                link = item['link']
                content = get_one_appmsgstat(link, headers, key, pass_ticket, appmsg_token)
                print(content)
                if 'appmsgstat' not in content:
                    sys.exit(1)
                read_num = str(content['appmsgstat'].get('read_num', 0))
                #read_num = str(content['appmsgstat'].get('real_read_num', 0))
                like_num = str(content['appmsgstat'].get('like_num', 0))
                #like_num = str(content['appmsgstat'].get('show_like', 0))
                old_like_num = str(content['appmsgstat'].get('old_like_num', 0))
                create_time = get_date(item['create_time'])
                # store information
                info = '"{}","{}","{}","{}","{}","{}","{}"'.format(aid,title,read_num,like_num,old_like_num,create_time,link)
                info_list.append(info)
    return info

if __name__ == '__main__':
    opts = get_opts()
    yaml_file =  opts.yaml
    json_file = opts.json

    config = read_config( yaml_file )
    #info = get_all_appmsgstat( json_file, config)
    info_list = get_all_article( json_file, config)
    with open("./authorst.txt", "w") as f:
        f.writelines("\n".join(info_list))
