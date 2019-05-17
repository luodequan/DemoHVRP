

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import nju.lzx.Constraint.*;
import nju.lzx.Constraint.MinimizeDistance.ConstraintData;
import nju.lzx.Data.*;
import nju.lzx.Algorithm.*;
import nju.lzx.Interface.*;
import nju.lzx.LocalSearchOperator.*;
import nju.lzx.InitialSolutionOperator.*;
import nju.lzx.VehicleReductionOperator.*;
import nju.lzx.Utility.*;
import nju.lzx.Route.*;


// TODO: Auto-generated Javadoc
/**
 * Heterogeneous Vehicle Routing Problem (HVRP)求解示例。
 */
public class HVRP {

	
	public static class InstanceHetero extends Instance{

		public int nt; //the number of types
		public int[] mt; //the number of available vehicles for each type
		public double[] Qt; //the capacity of each type of vehicle
		public double[] ut; //the unit travel cost of each type of vehicle
		public double[] ft; //the fixed cost of each type of vehicle
	}

	
	/**
	 * 主函数。
	 *
	 * @param args 参数。
	 * @throws IOException IO异常。
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		double t1 = System.nanoTime();
		InstanceHetero inst = (InstanceHetero) load_instance("data/Taillard_Heterogeneous/c100_20mix.txt");
		inst.m = 99;
		//设置模式参数
		inst.parameter.Mode.multi_thread_enable = true;
		//设置初始解参数
		inst.parameter.InitialSolution.log_print = false;
		//设置禁忌搜索参数
		inst.parameter.TabuSearch.maximum_iteration = 1000;
		inst.parameter.TabuSearch.maximum_tabu_tenure = 50;
		inst.parameter.TabuSearch.tenure_decay_rate = 0.99;
		inst.parameter.TabuSearch.mininum_tabu_tenure = 10;
		inst.parameter.TabuSearch.mininum_shake_tenure = 10;
		inst.parameter.TabuSearch.shake_tenure_increment = 10;
		inst.parameter.TabuSearch.maximum_shake_tenure = 500;
		inst.parameter.TabuSearch.minimum_shake_iteration = 100;
		inst.parameter.TabuSearch.log_print = true;
		inst.parameter.TabuSearch.log_detail = false;
		//设置算子参数
		inst.parameter.Operator.insertion_prune_threshhold = 1e6;
		inst.parameter.Operator.exchange_prune_threshhold = 1e6;
		inst.parameter.Operator.cross_prune_threshhold = 1e6;
		inst.parameter.Operator.remove_prune_threshhold = 0;
		inst.parameter.Operator.route_cross_threshhold = 1e6;
		
		
		//构造约束条件
		Constraint[] cnts = new Constraint[3];
		MinimizeDistance.ConstraintData[] dist_dats = new MinimizeDistance.ConstraintData[inst.nt];
		for(int i = 0; i < inst.nt; i++)
			dist_dats[i] = new MinimizeDistance.ConstraintData(inst.d, inst.ut[i]);
		cnts[0] = new MinimizeDistance(dist_dats, 0, 1);
		MinimizeFixedCost.ConstraintData[] fix_dats = new MinimizeFixedCost.ConstraintData[inst.nt];
		for(int i = 0; i < inst.nt; i++)
			fix_dats[i] = new MinimizeFixedCost.ConstraintData(inst.ft[i]);
		cnts[1] = new MinimizeFixedCost(fix_dats, 0, 1);
		CapacityConstraint.ConstraintData[] cap_dats = new CapacityConstraint.ConstraintData[inst.nt];
		for(int i = 0; i < inst.nt; i++)
			cap_dats[i] = new CapacityConstraint.ConstraintData(inst.q, inst.Qt[i]);
		cnts[2] = new CapacityConstraint(cap_dats, 0, true, 30);
		
		//构造算法算子
		Operator[] operators = new Operator[4];
		double[] coefs = new double[4];
		operators[0] = new RelocateBase(inst);
		coefs[0] = 1;
		operators[1] = new ExchangeBaseDeep(inst);
		coefs[1] = 0.5;
		operators[2] = new CrossBase(inst);
		coefs[2] = 1;
		operators[3] = new RelocateBaseIntra(inst);
		coefs[3] = 1;
		
		//构造需要访问的节点集合
		ArrayList<Atr> atrs = new ArrayList<Atr>();
		for(int i = 1; i < inst.n; i++){
			atrs.add(new Atr(i));
		}
		boolean[] exc = new boolean[inst.n];
		for(int i = 1; i < inst.n; i++)
			exc[i] = true;
		
		//构造初始解
		int[] tp = new int[inst.nt];
		int[] mt = new int[inst.nt];
		inst.mt[0] = 10000;
		inst.mt[1] = 5;
		inst.mt[2] = 0;
		for(int i = 0; i < inst.nt; i++){
			tp[i] = inst.nt - 1 - i;
			mt[i] = inst.mt[inst.nt - 1 - i];
		}
		GreedyGeneral greedy = new GreedyGeneral(inst, new GreedyGeneral.AlgorithmData(0, inst.nt, tp, mt), 0, 0, new InsertBase(inst, cnts));
		ArrayList<Route> s = greedy.generate(atrs);
		System.out.println(greedy.toString(s));
		greedy.check(s, true, exc);
		System.out.println("feasibility of the initial solution>>>" + greedy.is_feasible(s) + "\t" + s.size() + "\t" + greedy.get_total_cost(s));
		//System.exit(0);
		
		//最小化行驶距离
		TabuSearch tabu = new TabuSearch(inst, operators, coefs);
		s = toDeep(inst, s);
		tabu.check(s, true, exc);
		for(int i = 0; i < s.size(); i++){
			s.get(i).relax(2);
		}
		ArrayList<Route> bs = tabu.solve(s);

		tabu.check(bs, true, exc);
		//System.out.println(tabu.toString(bs));
		double t2 = System.nanoTime();
		System.out.println("computation time>>>" + (t2 - t1) / 1e9);
		//inst.statistics.print();
	}
	
	/**
	 * 加载算例。
	 *
	 * @param path 算例路径。
	 * @return 返回算例。
	 * @throws FileNotFoundException IO异常。
	 */
	public static Instance load_instance(String path) throws FileNotFoundException{
		InstanceHetero inst = new InstanceHetero();
		Scanner cin = new Scanner(new BufferedReader(new FileReader(path)));
		inst.n = cin.nextInt() + 1;
		inst.q = new double[inst.n];
		inst.s = new double[inst.n];
		inst.lng = new double[inst.n];
		inst.lat = new double[inst.n];
		for(int i = 0; i < inst.n; i++){
			cin.next();
			inst.lng[i] = cin.nextDouble();
			inst.lat[i] = cin.nextDouble();
			inst.q[i] = cin.nextDouble();
		}
		String str = cin.nextLine();
		for(int i = 0; i < 4; i++)
			cin.nextLine();
		ArrayList<Double> Q_array = new ArrayList<Double>();
		ArrayList<Double> f_array = new ArrayList<Double>();
		ArrayList<Double> u_array = new ArrayList<Double>();
		ArrayList<Integer> m_array = new ArrayList<Integer>();
		while(cin.hasNextLine()){
			String line = cin.nextLine();
			if(line.length() < 2)
				break;
			String[] sub = line.split(" ");
			Q_array.add(Double.parseDouble(sub[0]));
			f_array.add(Double.parseDouble(sub[1]));
			u_array.add(Double.parseDouble(sub[2]));
			m_array.add(Integer.parseInt(sub[3]));
		}
		cin.close();
		inst.nt = Q_array.size();
		inst.Qt = new double[inst.nt];
		inst.ut = new double[inst.nt];
		inst.ft = new double[inst.nt];
		inst.mt = new int[inst.nt];
		inst.d = new double[inst.n][inst.n];
		inst.t = new double[inst.n][inst.n];
		for(int i = 0; i < inst.nt; i++){
			inst.Qt[i] = Q_array.get(i);
			inst.ft[i] = f_array.get(i);
			//inst.ut[i] = u_array.get(i);
			inst.ut[i] = 1;
			//inst.mt[i] = m_array.get(i);
			inst.mt[i] = 100000;
		}
		inst.mt[1] = 5;
		inst.mt[2] = 0;
		for(int i = 0; i < inst.n; i++){
			for(int j = i + 1; j < inst.n; j++){
				inst.d[i][j] = inst.d[j][i] = inst.t[i][j] = inst.t[j][i] = 
						Math.sqrt((inst.lng[i] - inst.lng[j]) * (inst.lng[i] - inst.lng[j]) + (inst.lat[i] - inst.lat[j]) * (inst.lat[i] - inst.lat[j]));
			}
		}
		return inst;
	}

	
	/**
	 * 将一个普通解转换成一个可以被ExchangeBaseDeep和RelocateBseIntra操作的解，即路径中包含子路径。
	 *
	 * @param inst 算例信息。
	 * @param s 当前解。
	 * @return 返回新的解。
	 */
	public static ArrayList<Route> toDeep(Instance inst, ArrayList<Route> s){
		ArrayList<Route> ns = new ArrayList<Route>();
		for(int i = 0; i < s.size(); i++){
			Route r = s.get(i);
			Reference ref = r.get_reference();
			Route nr = new RouteBase(inst, ref.len, ref.seq, false, true);
			Constraint[] _cnts = r.get_constraints();
			Constraint[] _cnts2 = new Constraint[_cnts.length];
			for(int j = 0; j < _cnts.length; j++){
				_cnts2[j] = _cnts[j].copy(nr.get_reference());
			}
			nr.add_constraints(_cnts2);
			ns.add(nr);
		}
		return ns;
	}
}
