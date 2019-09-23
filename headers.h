// -----------------------------------------------------------------------------
//  Copyright (c) 2013-2018 Sun Yat-Sen University (SYSU). All Rights Reserved.
//
//  SYSU grants permission to use, copy, modify, and distribute this software
//  and its documentation for NON-COMMERCIAL purposes and without fee, provided 
//  that this copyright notice appears in all copies.
//
//  SYSU provides this software "as is," without representations or warranties
//  of any kind, either expressed or implied, including but not limited to the
//  implied warranties of merchantability, fitness for a particular purpose, 
//  and noninfringement. SYSU shall not be liable for any damages arising from
//  any use of this software.
//
//  Author: Qiang Huang (huangq2011@gmail.com)
//
//  Release Date:  23-04-2018
//  Last Modified: 23-04-2018
//  Version 1.0.0
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <errno.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"
#include "qdafn.h"
#include "drusilla_select.h"
#include "rqalsh.h"
#include "rqalsh_star.h"
#include "ml_rqalsh.h"
#include "afn.h"

using namespace std;

